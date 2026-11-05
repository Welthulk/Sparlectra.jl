# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# file: src/stateestimation/measurement_placement.jl
# purpose: measurement and PMU placement (concept 1.1.0, section 1.3):
#          suggest_measurement_placement proposes meter locations step by
#          step until a set is observable or has no critical measurement,
#          on the decoupled P-theta and Q-V models of the observability
#          check, per network island; print_measurement_placement prints
#          the step list.

# the values of the `kind` and `target` keywords
const _PLACEMENT_KINDS = (:injection, :flow, :pmu)
const _PLACEMENT_TARGETS = (:observable, :no_critical)

## The candidates of one network island (its net `inet`, the active rows
## `meas` in its numbering, the working net `prep` and its active
## non-current rows `act`) for `kind`. Each candidate is
## `(quantity, bus, branch, side, rows)`: `bus` in the numbering of `inet`
## (0 for a flow), `branch` the branch of the working net (0 otherwise),
## `side` the measured end, `rows` the rows it adds in the numbering of the
## working net.
## - `:injection`: P, Q and U (a Vm row) at a bus without an active row of
##   that type, each one candidate with one row (P closes the P-theta model,
##   Q and U the Q-V model).
## - `:pmu`: a `VaMeas` at a bus without one.
## - `:flow`: the P and the Q flow at one end of a closed branch that has
##   neither an active P nor an active Q flow row at that end (a meter at
##   an end measures both), one candidate with two rows.
## A bus inside a closed-link cluster of more than one bus gets no
## candidate (its rows are the cluster's, as in the restoration), nor does
## a bus the contraction drops. Rows carry value 0 (1 pu for U) and the
## default sigma: only their place enters the observability check.
function _placement_candidates(inet::Net, meas::Vector{Measurement}, prep::SEWorkingNet, act::Vector{Measurement}, kind::Symbol)
  sd = measurementStdDevs()
  snet = prep.snet
  cands = NamedTuple[]
  if kind === :flow
    has = Set{Tuple{Int,Symbol}}()
    for m in act
      (m.linkIdx === nothing && m.branchIdx !== nothing && m.typ in (PflowMeas, QflowMeas)) || continue
      push!(has, (m.branchIdx, m.direction == :to ? :to : :from))
    end
    for (k, br) in enumerate(snet.branchVec)
      _branch_terminal_state(br) == :closed || continue
      Int(br.fromBus) == Int(br.toBus) && continue
      for side in (:from, :to)
        (k, side) in has && continue
        rows = [Measurement(typ = PflowMeas, value = 0.0, sigma = sd[PflowMeas], branchIdx = k, direction = side, id = "placement_candidate"), Measurement(typ = QflowMeas, value = 0.0, sigma = sd[QflowMeas], branchIdx = k, direction = side, id = "placement_candidate")]
        push!(cands, (quantity = :PQ, bus = 0, branch = k, side = side, rows = rows))
      end
    end
    return cands
  end
  nb = length(inet.nodeVec)
  reps = _active_link_representative_map(inet)
  clusterSize = Dict{Int,Int}()
  for b = 1:nb
    clusterSize[reps[b]] = get(clusterSize, reps[b], 0) + 1
  end
  hasRow = Set{Tuple{MeasurementType,Int}}()
  for m in meas
    (m.active && m.busIdx !== nothing && m.linkIdx === nothing) && push!(hasRow, (m.typ, m.busIdx))
  end
  types = kind === :pmu ? ((:Va, VaMeas),) : ((:P, PinjMeas), (:Q, QinjMeas), (:U, VmMeas))
  for b = 1:nb
    clusterSize[reps[b]] > 1 && continue
    s = get(prep.busmap, prep.reps[b], 0)
    s == 0 && continue
    for (q, typ) in types
      (typ, b) in hasRow && continue
      push!(cands, (quantity = q, bus = b, branch = 0, side = :none, rows = [Measurement(typ = typ, value = typ == VmMeas ? 1.0 : 0.0, sigma = sd[typ], busIdx = s, id = "placement_candidate")]))
    end
  end
  return cands
end

## The state of one network island during the placement: the models of
## the active rows `act` on the working net, their null spaces and the
## observable islands. `aug` is the P-theta model in the columns the
## candidates are tested in: with the PMU offset column whenever a PMU
## candidate would carry it (`:pmu` under `pmu_ref_offset = auto`), so the
## first PMU of a set without one is tested against an empty offset column
## (a null direction of its own); the deficits reported are those of the
## set itself (`angle_d` with the offset only when the set has a PMU).
## `unobs` are the unobservable branches of the full boundary iteration in
## the numbering of `inet`.
function _placement_state(inet::Net, prep::SEWorkingNet, act::Vector{Measurement}, kind::Symbol, mode::Symbol, dense_max_n::Integer)
  snet = prep.snet
  wReal = _va_offset_active(act, mode)
  wAug = kind === :pmu ? mode === :auto : wReal
  pt = _ptheta_model(snet, act, wReal)
  pn = _decoupled_null(pt.A, dense_max_n)
  pnAug = wAug == wReal ? pn : _decoupled_null(_ptheta_model(snet, act, wAug).A, dense_max_n)
  qt = _qv_model(snet, act)
  qn = _decoupled_null(qt.A, dense_max_n)
  piece = _observable_islands_working(inet, prep, act, pt, pn, qn, dense_max_n, typemax(Int))
  return (pnAug = pnAug, qn = qn, angle_d = pn.d, mag_d = qn.d, unobs = sort!(unique(piece.unobs)), slackIdx = pt.slackIdx)
end

## The number of branches of `edges` whose ends lie in different groups of
## one plain pass over the null basis `N` of the P-theta model (columns
## normalized by `cn`, probe tolerance `tau`; `_angle_island_labels`, no
## boundary iteration): the ordering key of the placement, cheap enough to
## evaluate per candidate. Bus rows are the angle columns; the circulation
## and offset columns are no bus.
function _placement_plain_unobservable(N::Matrix{Float64}, cn::Vector{Float64}, tau::Float64, nb::Int, slackIdx::Int, edges::Vector{Tuple{Int,Int,Int}})::Int
  Wb = zeros(Float64, nb, size(N, 2))
  cnb = fill(Inf, nb)
  for i = 1:nb
    i == slackIdx && continue
    p = i < slackIdx ? i : i - 1
    Wb[i, :] = view(N, p, :)
    cnb[i] = cn[p]
  end
  labels = _angle_island_labels(Wb, cnb, edges, tau)
  return count(e -> labels[e[2]] != labels[e[3]], edges)
end

## The critical rows (numerical and structural, positions in `rows`) of the
## rows `rows` of the AC Jacobian `Hall` (all rows of the network island and
## of its candidates), judged by `_observability_core` as the verdict judges
## them: rows then columns normalized, the FD-aware rank tolerance, the
## configured criticality and rank methods. `offsetCol` is the PMU offset
## column of `Hall` (0 without), dropped when `hasVa` is false (the set
## carries no offset then). -1 when the rows are not observable in the AC
## Jacobian.
function _placement_critical(Hall::SparseMatrixCSC{Float64}, rows::Vector{Int}, offsetCol::Int, hasVa::Bool, cfg::StateEstimationConfig, dense_max_n::Integer)::Int
  cols = (offsetCol == 0 || hasVa) ? collect(1:size(Hall, 2)) : [j for j = 1:size(Hall, 2) if j != offsetCol]
  Hs, _ = _row_column_normalized(Hall[rows, cols])
  tol = cfg.rank_tol_factor * cfg.jac_eps * _sigma_max(Hs)
  core = _observability_core(Hs, collect(1:length(rows)); tol = tol, criticality_method = cfg.criticality_method, rank_method = cfg.rank_method, dense_max_n = dense_max_n)
  core.numerical_observable || return -1
  return length(union(core.numerical_critical_measurement_indices, core.structural_critical_measurement_indices))
end

"""
    suggest_measurement_placement(net, measurements; kind, target = :observable, max_steps = state_estimation.observability.placement_max_steps, dense_max_n = $(_SE_DENSE_LINALG_MAX_N)) -> NamedTuple
    suggest_measurement_placement(net; kwargs...) -> NamedTuple

Where to place the next meters so that the set becomes observable
(`target = :observable`) or keeps no critical measurement
(`target = :no_critical`). One meter per step, per network island, never
across a network island border; a network island without measurements is
left out. The caller's vector is not changed.

Candidates by `kind` (one candidate per place and quantity):
- `:injection`: P, Q or U (a Vm row) at a bus without an active row of
  that type, one row each. P closes the P-theta model, Q and U the Q-V
  model.
- `:flow`: P and Q at one end of a closed branch without an active P or Q
  flow row at that end, two rows.
- `:pmu`: a voltage angle (`VaMeas`) at a bus without one, one row. Under
  `pmu_ref_offset = auto` the PMU time base is an extra state, so the
  first PMU of a set without one fixes only that offset and lowers no
  deficit (a PMU at the slack fixes nothing else either); the placement
  counts that step.
No candidate on a bus inside a closed-link cluster of more than one bus.

Target `:observable`: steps until the P-theta and the Q-V deficit of
[`evaluate_global_observability`](@ref) are both 0. A candidate is valid
when it lowers the deficit of its model by one (the rank test of
[`restore_observability`](@ref), on the null space of the set plus the
rows taken so far). With one row per candidate (`:injection`, `:pmu`)
every valid sequence has as many steps as the deficits sum to (plus the
offset step of `:pmu`), so the number is the minimum; the order is the
choice: the valid candidate after which one plain null-space pass finds
the fewest unobservable branches (the most removed), ties by bus, then P,
Q, U. With `:flow` the candidate with the largest rank gain (P-theta plus
Q-V) is taken, ties by the same branch count, then branch and end; there
is no claim of a minimum. A magnitude gap no candidate of the kind closes
(`:pmu`, or `:flow` without a magnitude anchor) stops with
`stop_reason = :no_candidate`; under `:pmu` without a step when the
P-theta deficit is 0.

Target `:no_critical`: first as `:observable`, then steps until no row of
the network island is critical in the AC Jacobian (numerical or
structural, as `numerical_critical_measurement_indices` and
`structural_critical_measurement_indices` of
[`evaluate_global_observability`](@ref)): the candidate after which the
fewest rows are critical, ties by bus (branch and end), then quantity; a
step must lower that number (the offset step of `:pmu` excepted), else
`:no_candidate`; a set observable in the decoupled models but not in the
AC Jacobian stops with `:ac_rank_deficit`. Each trial is one
criticality computation, so the cost grows with candidates times steps.

`max_steps` caps the number of steps over all network islands (`0`
allowed); reaching it before the target sets `stop_reason = :max_steps`.
With several network islands `stop_reason` is the one of the first that
does not reach the target.
An abort of the surrounding job is checked between steps and every 64
candidates.

Returns `(kind, target, max_steps, reached, stop_reason, steps, rows,
network_islands, angle_rank_deficit_before, magnitude_rank_deficit_before,
angle_rank_deficit_after, magnitude_rank_deficit_after,
unobservable_branches_before, unobservable_branches_after, seconds)`:
- `steps`: one row per step `(step, network_island, phase, kind, quantity,
  bus, bus_name, branch, side, location, angle_rank_deficit_before,
  angle_rank_deficit_after, magnitude_rank_deficit_before,
  magnitude_rank_deficit_after, unobservable_branches_before,
  unobservable_branches_after, critical_before, critical_after, seconds)`.
  `phase` is `:observable` or `:no_critical`; `quantity` `:P`, `:Q`, `:U`,
  `:Va` or `:PQ` (a flow); `bus` (0 for a flow) and `branch` (0 otherwise)
  in the caller's numbering, `side` `:from`, `:to` or `:none`; the
  unobservable branches are counted by the boundary iteration of
  [`identify_observable_islands`](@ref) before and after the step; the
  critical counts are `nothing` in the phase `:observable`; `seconds` is
  the run time of the step.
- `rows`: the suggested `Measurement`s in the caller's numbering (id
  `placement_<default id>`, value 0, 1 pu for U, the default sigma of
  `measurementStdDevs`): the places only, for an observability check, no
  readings.
- `network_islands`: one row per measured network island
  `(network_island, n_steps, angle_rank_deficit_before,
  magnitude_rank_deficit_before, angle_rank_deficit_after,
  magnitude_rank_deficit_after, reached)`.
- the summed deficits and the unobservable branches (caller's numbering)
  before and after, and the total run time.

The step times are in the result; the SE service writes them to
performance.log. Errors like `evaluate_global_observability` when no
active measurement is left, and on an unknown `kind` or `target` or a
negative `max_steps`.
"""
function suggest_measurement_placement(net::Net, measurements::Vector{Measurement}; kind::Symbol, target::Symbol = :observable, max_steps::Integer = state_estimation_config().observability.placement_max_steps, dense_max_n::Integer = _SE_DENSE_LINALG_MAX_N)
  kind in _PLACEMENT_KINDS || error("suggest_measurement_placement: kind must be one of $(join(_PLACEMENT_KINDS, ", ")) (got $(kind))")
  target in _PLACEMENT_TARGETS || error("suggest_measurement_placement: target must be one of $(join(_PLACEMENT_TARGETS, ", ")) (got $(target))")
  max_steps >= 0 || error("suggest_measurement_placement: max_steps must not be negative (got $(max_steps))")
  t0 = time_ns()
  cfg = state_estimation_config()
  mode = cfg.pmu_ref_offset
  part = _se_island_partition(net, measurements)
  names = _bus_name_by_idx(net)
  steps = NamedTuple[]
  rows = Measurement[]
  perIsland = NamedTuple[]
  unobsBefore = Int[]
  unobsAfter = Int[]
  stop = :reached
  for (k, nrow) in enumerate(part.rows)
    isempty(part.perIsland[k]) && continue
    # the network island on its own net, as in restore_observability
    inet, meas_k, toCallerBus, toCallerBranch = if length(part.rows) == 1
      (net, measurements, collect(1:length(net.nodeVec)), collect(1:length(net.branchVec)))
    else
      sub = _se_island_subnet(net, nrow)
      m_k, _ = _se_remap_island_measurements(sub, measurements, part.perIsland[k])
      (sub.inet, m_k, collect(nrow.buses), collect(nrow.branches))
    end
    # a network island whose rows the working net does not keep (no
    # active non-current row left) is left out like an unmeasured one
    prep = _se_prepare(inet, meas_k)
    act0 = Measurement[m for m in prep.meas if !(m.typ in (ImagMeas, IaMeas))]
    isempty(act0) && continue
    snet = prep.snet
    nbw = length(snet.nodeVec)
    act = copy(act0)
    edges = Tuple{Int,Int,Int}[]
    for (j, br) in enumerate(snet.branchVec)
      _branch_terminal_state(br) == :closed || continue
      f = Int(br.fromBus)
      t = Int(br.toBus)
      f == t || push!(edges, (j, f, t))
    end
    cands = _placement_candidates(inet, meas_k, prep, act, kind)
    # the candidate rows in both models, once: flat list, owner per row;
    # the columns do not change during the steps (the offset column of
    # `aug` is fixed for the kind)
    flat = Measurement[]
    owner = Int[]
    for (i, c) in enumerate(cands)
      for r in c.rows
        push!(flat, r)
        push!(owner, i)
      end
    end
    wAug = kind === :pmu ? mode === :auto : _va_offset_active(act, mode)
    Mp = _decoupled_model(snet, flat, :ptheta, wAug)
    Mq = _decoupled_model(snet, flat, :qv, false)
    MpT = sparse(transpose(Mp.A))
    MqT = sparse(transpose(Mq.A))
    pRow = zeros(Int, length(cands))
    qRow = zeros(Int, length(cands))
    for (r, q) in enumerate(Mp.rows)
      pRow[owner[q]] = r
    end
    for (r, q) in enumerate(Mq.rows)
      qRow[owner[q]] = r
    end
    open = trues(length(cands))
    st = _placement_state(inet, prep, act, kind, mode, dense_max_n)
    first0 = st
    append!(unobsBefore, Int[toCallerBranch[j] for j in st.unobs])
    nSteps = 0
    reachedK = false
    # why this network island stopped; the result keeps the first reason
    stopK = :reached
    # the caller's view of one candidate
    function place(c)
      if c.branch == 0
        b = toCallerBus[c.bus]
        return (bus = b, bus_name = get(names, b, string(b)), branch = 0, location = get(names, b, string(b)))
      end
      kc = toCallerBranch[prep.branches[c.branch]]
      br = net.branchVec[kc]
      return (bus = 0, bus_name = "", branch = kc, location = string(get(names, Int(br.fromBus), string(br.fromBus)), " - ", get(names, Int(br.toBus), string(br.toBus)), " (", c.side, " end)"))
    end
    # takes candidate i: its rows join the set (working and caller
    # numbering), the state is computed again in full
    function takestep!(i, phase, t1, critBefore, critAfter)
      c = cands[i]
      open[i] = false
      append!(act, c.rows)
      p = place(c)
      for r in c.rows
        if c.branch == 0
          push!(rows, Measurement(typ = r.typ, value = r.value, sigma = r.sigma, busIdx = p.bus, id = string("placement_", _default_measurement_id(r.typ; busIdx = p.bus))))
        else
          push!(rows, Measurement(typ = r.typ, value = r.value, sigma = r.sigma, branchIdx = p.branch, direction = c.side, id = string("placement_", _default_measurement_id(r.typ; branchIdx = p.branch, direction = c.side))))
        end
      end
      before = st
      st = _placement_state(inet, prep, act, kind, mode, dense_max_n)
      nSteps += 1
      push!(steps, (step = length(steps) + 1, network_island = k, phase = phase, kind = kind, quantity = c.quantity, bus = p.bus, bus_name = p.bus_name, branch = p.branch, side = c.side, location = p.location,
        angle_rank_deficit_before = before.angle_d, angle_rank_deficit_after = st.angle_d, magnitude_rank_deficit_before = before.mag_d, magnitude_rank_deficit_after = st.mag_d,
        unobservable_branches_before = length(before.unobs), unobservable_branches_after = length(st.unobs), critical_before = critBefore, critical_after = critAfter, seconds = (time_ns() - t1) / 1.0e9))
      return nothing
    end
    qtyRank(q) = q === :P ? 1 : (q === :Q ? 2 : (q === :U ? 3 : 4))
    # phase 1: observable
    while true
      # the deficits of the set itself decide: `pnAug` of a set without a
      # PMU also counts the empty offset column of the PMU candidates
      if st.angle_d + st.mag_d == 0
        reachedK = true
        break
      end
      length(steps) >= max_steps && (stopK = :max_steps; break)
      # a PMU only enters the P-theta model: with its deficit closed it
      # cannot lower a magnitude gap, and the first PMU of a set without one
      # would only fix the offset
      (kind === :pmu && st.angle_d == 0) && (stopK = :no_candidate; break)
      sparlectra_check_abort()
      t1 = time_ns()
      pA = st.pnAug
      qA = st.qn
      plain0 = _placement_plain_unobservable(pA.N, pA.cn, pA.tau, nbw, st.slackIdx, edges)
      best = 0
      bestKey = (0, 0, 0, 0)
      for (i, c) in enumerate(cands)
        open[i] || continue
        i % 64 == 0 && sparlectra_check_abort()
        gp, cp = pRow[i] == 0 ? (false, Float64[]) : _null_row_share(pA.N, pA.cn, pA.tau, Vector(MpT[:, pRow[i]]))
        gq, _ = qRow[i] == 0 ? (false, Float64[]) : _null_row_share(qA.N, qA.cn, qA.tau, Vector(MqT[:, qRow[i]]))
        gain = Int(gp) + Int(gq)
        gain == 0 && continue
        # the null space with the row added (as in `_restore_select`)
        plain = gp ? _placement_plain_unobservable(pA.N * nullspace(reshape(cp, 1, :)), pA.cn, pA.tau, nbw, st.slackIdx, edges) : plain0
        # smaller is better: largest gain (flow only), fewest unobservable
        # branches after one plain pass, then the place
        key = c.branch == 0 ? (0, plain, c.bus, qtyRank(c.quantity)) : (-gain, plain, c.branch, c.side === :from ? 1 : 2)
        if best == 0 || key < bestKey
          best = i
          bestKey = key
        end
      end
      best == 0 && (stopK = :no_candidate; break)
      takestep!(best, :observable, t1, nothing, nothing)
    end
    # phase 2: no critical row
    if reachedK && target === :no_critical
      reachedK = false
      x = _initial_state_vector(snet, st.slackIdx; flatstart = cfg.flatstart, withVaOffset = _va_offset_active(vcat(act, flat), mode), rotate_behind_pst = cfg.start_mode.rotate_behind_pst)
      allRows = vcat(act, flat)
      wAll = _va_offset_active(allRows, mode)
      Hall, _ = _measurement_jacobian_fd(allRows, snet, x, st.slackIdx, nbw, createYBUS(net = snet); eps = cfg.jac_eps, withVaOffset = wAll)
      size(Hall, 1) == length(allRows) || error("suggest_measurement_placement: the AC Jacobian has $(size(Hall, 1)) rows for $(length(allRows)) measurement rows")
      Hall = sparse(Hall)
      offsetCol = wAll ? size(Hall, 2) : 0
      # rows of Hall per candidate: behind the rows of the set at the start
      # of this phase
      nAct = length(act)
      candRows = [Int[] for _ in cands]
      for (q, i) in enumerate(owner)
        push!(candRows[i], nAct + q)
      end
      cur = collect(1:nAct)
      curVa = any(m -> m.typ == VaMeas, act)
      crit = _placement_critical(Hall, cur, offsetCol, curVa, cfg, dense_max_n)
      # observable in the decoupled models but not in the AC Jacobian: no
      # criticality to lower
      crit < 0 && (stopK = :ac_rank_deficit)
      while crit >= 0
        if crit == 0
          reachedK = true
          break
        end
        length(steps) >= max_steps && (stopK = :max_steps; break)
        sparlectra_check_abort()
        t1 = time_ns()
        best = 0
        # the first PMU of a set without one only fixes the offset (its own
        # row is critical), so that step is taken without lowering the
        # count, as in the first phase
        offsetStep = kind === :pmu && offsetCol != 0 && !curVa
        bestCrit = offsetStep ? typemax(Int) : crit
        for (i, c) in enumerate(cands)
          open[i] || continue
          i % 64 == 0 && sparlectra_check_abort()
          cv = _placement_critical(Hall, vcat(cur, candRows[i]), offsetCol, curVa || kind === :pmu, cfg, dense_max_n)
          # a trial must lower the count; candidates come by bus (branch
          # and end), then quantity, so the first with the fewest wins
          (cv >= 0 && cv < bestCrit) || continue
          best = i
          bestCrit = cv
        end
        best == 0 && (stopK = :no_candidate; break)
        append!(cur, candRows[best])
        curVa = curVa || kind === :pmu
        takestep!(best, :no_critical, t1, crit, bestCrit)
        crit = bestCrit
      end
    end
    (stop === :reached && !reachedK) && (stop = stopK)
    append!(unobsAfter, Int[toCallerBranch[j] for j in st.unobs])
    push!(perIsland, (network_island = k, n_steps = nSteps, angle_rank_deficit_before = first0.angle_d, magnitude_rank_deficit_before = first0.mag_d, angle_rank_deficit_after = st.angle_d, magnitude_rank_deficit_after = st.mag_d, reached = reachedK))
  end
  isempty(perIsland) && error("suggest_measurement_placement: no active measurements (current magnitudes are excluded; the system must be observable without them)")
  reached = all(r -> r.reached, perIsland)
  reached && (stop = :reached)
  return (kind = kind, target = target, max_steps = Int(max_steps), reached = reached, stop_reason = stop, steps = steps, rows = rows, network_islands = perIsland,
    angle_rank_deficit_before = sum(r -> r.angle_rank_deficit_before, perIsland), magnitude_rank_deficit_before = sum(r -> r.magnitude_rank_deficit_before, perIsland),
    angle_rank_deficit_after = sum(r -> r.angle_rank_deficit_after, perIsland), magnitude_rank_deficit_after = sum(r -> r.magnitude_rank_deficit_after, perIsland),
    unobservable_branches_before = sort!(unique(unobsBefore)), unobservable_branches_after = sort!(unique(unobsAfter)), seconds = (time_ns() - t0) / 1.0e9)
end

suggest_measurement_placement(net::Net; kwargs...) = suggest_measurement_placement(net, Measurement[m for m in net.measurements]; kwargs...)

"""
    print_measurement_placement([io,] result)

Print the result of [`suggest_measurement_placement`](@ref): one header
line (kind, target, number of steps, whether the target was reached and
why it stopped), then one line per step with its network island, the
place and quantity, the angle and magnitude deficits and the unobservable
branches before and after, the critical rows in the phase `:no_critical`,
and the run time.
"""
function print_measurement_placement(io::IO, r::NamedTuple)
  n = length(r.steps)
  status = r.reached ? "target reached" : r.stop_reason === :max_steps ? string("stopped at max_steps = ", r.max_steps, ", target not reached") :
           r.stop_reason === :ac_rank_deficit ? "stopped, the AC Jacobian has a rank deficit the decoupled models do not show" : "stopped, no candidate of this kind lowers the gap"
  println(io, "measurement placement: kind ", r.kind, ", target ", r.target, ", ", n, n == 1 ? " step" : " steps", ", ", status, "; deficit angles ", r.angle_rank_deficit_before, " -> ", r.angle_rank_deficit_after, ", magnitudes ", r.magnitude_rank_deficit_before, " -> ", r.magnitude_rank_deficit_after)
  what(s) = s.quantity === :PQ ? string("P and Q flow at ", s.location) : string(s.quantity === :Va ? "PMU angle" : (s.quantity === :U ? "U" : string(s.quantity, " injection")), " at bus ", s.location)
  for s in r.steps
    crit = s.critical_before === nothing ? "" : string(", critical rows ", s.critical_before, " -> ", s.critical_after)
    println(io, "  step ", s.step, " (network island ", s.network_island, "): ", what(s), ", deficit angles ", s.angle_rank_deficit_before, " -> ", s.angle_rank_deficit_after, ", magnitudes ", s.magnitude_rank_deficit_before, " -> ", s.magnitude_rank_deficit_after,
      ", unobservable branches ", s.unobservable_branches_before, " -> ", s.unobservable_branches_after, crit, ", ", round(s.seconds; sigdigits = 3), " s")
  end
  return nothing
end

print_measurement_placement(r::NamedTuple) = print_measurement_placement(stdout, r)
