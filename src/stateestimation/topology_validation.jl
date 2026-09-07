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

# file: src/stateestimation/topology_validation.jl
# purpose: topology validation for the state estimation (0.10.0):
#          stage 1 pre-checks (validate_topology, pure and linear), and
#          stage 3 hypothesis testing on WORKING COPIES
#          (test_topology_hypotheses). Stage 2, the post-SE fingerprint
#          classification, lives in runse_diagnostics. All three stages are
#          ADVISORY: they produce findings and warnings, they never mutate
#          a status, a measurement, or the model.

## station of a bus: the representative of its closed-link contraction
## cluster (a bus without links is its own station). Findings and the
## stage-2 clustering aggregate over stations, not raw buses.
function _topology_station_map(net::Net)
  reps = _active_link_representative_map(net)
  members = Dict{Int,Vector{Int}}()
  for b in eachindex(reps)
    push!(get!(members, reps[b], Int[]), b)
  end
  return reps, members
end

## human-readable location label of a station
function _topology_station_label(net::Net, rep::Int, members::Dict{Int,Vector{Int}})
  nby = _bus_name_by_idx(net)
  ms = get(members, rep, [rep])
  base = get(nby, rep, string(rep))
  return length(ms) > 1 ? string(base, " (+", length(ms) - 1, " linked)") : base
end

## severity from how far the evidence exceeds its threshold (open design
## decision 1): below 2x the k factor :warning, at or above :strong. The
## dead-branch check never exceeds :warning (legitimately unloaded branches
## exist; the check only fires with loaded neighbours, see below).
_topology_severity(ratio::Float64, k::Float64) = ratio >= 2.0 * k ? :strong : :warning

"""
    validate_topology(net, measurements = net.measurements; kwargs...) -> NamedTuple

Stage-1 topology pre-checks: pure, linear in the measurement
count, no state estimation involved. Returns
`(findings, summary, n_checked_branches, n_checked_links)` where each
finding is `(stage = :precheck, kind, location, evidence, severity)`.
ADVISORY only: nothing is blocked or mutated; the caller decides what to
do with the findings. Out-of-service elements are exempt from the
plausibility checks; only the status-contradiction check
(`:open_element_with_flow`) looks at them, because a measured flow over an
open element IS the topology error.

Checks (thresholds as sigma multiples, keywords mirror the
`state_estimation.topology_*` configuration):

- `:open_element_with_flow` (`k_open`, default 4.0): an OPEN branch or
  link carries an active flow/current measurement with `|value| > k * sigma`.
- `:closed_element_without_flow` (`k_dead`, default 3.0): a CLOSED branch
  whose present flow/current measurements are all below `k * sigma` at
  both ends, while at least one terminal station carries other measured
  flows above the threshold (lower severity; a legitimately unloaded
  branch with equally quiet neighbours never fires).
- `:closed_link_voltage_mismatch` (`k_v`, default 4.0): a closed link with
  voltage-magnitude measurements on both sides disagreeing by more than
  `k * sqrt(sigma_a^2 + sigma_b^2)`.
- `:kcl_violation` (`k_kcl`, default 4.0): a COMPLETELY measured node
  (injection measurement present, every closed attached branch carries a
  flow measurement at this end) whose balance misses by more than
  `k * sqrt(sum sigma^2)`, for P and Q separately. Partially measured
  nodes and nodes inside closed-link clusters are skipped, never guessed.
"""
function validate_topology(net::Net, measurements::Vector{Measurement} = Measurement[m for m in net.measurements]; k_open::Float64 = state_estimation_config().topology_open_flow_k, k_dead::Float64 = state_estimation_config().topology_dead_flow_k, k_v::Float64 = state_estimation_config().topology_voltage_k, k_kcl::Float64 = state_estimation_config().topology_kcl_k)
  findings = NamedTuple[]
  active = [m for m in measurements if m.active]
  nby = _bus_name_by_idx(net)
  busname(b) = get(nby, Int(b), string(Int(b)))
  reps, members = _topology_station_map(net)

  flowTypes = (PflowMeas, QflowMeas, ImagMeas)
  isflow(m) = m.typ in flowTypes && m.branchIdx !== nothing
  # end bus a branch-flow row measures at (direction :from / :to)
  endbus(m) = begin
    br = net.branchVec[m.branchIdx]
    m.direction == :to ? Int(br.toBus) : Int(br.fromBus)
  end

  # Row indices, built once in O(rows). Every check below asks either "the
  # rows of element k" or "the row of type t measured at (branch k, bus b)".
  # Answering that by scanning `active` per element makes each check
  # O(elements x rows): on a 25000-bus network with 75000 rows that is
  # billions of comparisons and was measured as 45 s of a 55 s estimation
  # (the precheck, not the estimator, dominated the run). The maps below
  # preserve the original iteration order, so the findings are unchanged.
  flows_by_branch = Dict{Int,Vector{Measurement}}()
  for m in active
    isflow(m) || continue
    push!(get!(Vector{Measurement}, flows_by_branch, Int(m.branchIdx)), m)
  end
  flows_by_link = Dict{Int,Vector{Measurement}}()
  for m in active
    (m.linkIdx !== nothing && m.typ in flowTypes) || continue
    push!(get!(Vector{Measurement}, flows_by_link, Int(m.linkIdx)), m)
  end
  # closed branches attached to a bus (a self-loop is listed once)
  closed_branches_by_bus = Dict{Int,Vector{Int}}()
  for (k, br) in enumerate(net.branchVec)
    _branch_terminal_state(br) == :closed || continue
    fb = Int(br.fromBus)
    tb = Int(br.toBus)
    push!(get!(Vector{Int}, closed_branches_by_bus, fb), k)
    fb == tb || push!(get!(Vector{Int}, closed_branches_by_bus, tb), k)
  end
  # (branch, measured end bus, type) -> FIRST matching position in `active`,
  # which is exactly what the findfirst in check 4 returned
  first_flow_row = Dict{Tuple{Int,Int,MeasurementType},Int}()
  for (i, m) in enumerate(active)
    (m.branchIdx !== nothing && m.typ in flowTypes) || continue
    key = (Int(m.branchIdx), endbus(m), m.typ)
    haskey(first_flow_row, key) || (first_flow_row[key] = i)
  end
  # in-service shunts per bus (check 4 reads them per injection row)
  shunts_by_bus = Dict{Int,Vector{eltype(net.shuntVec)}}()
  for sh in net.shuntVec
    sh.status == 1 || continue
    push!(get!(Vector{eltype(net.shuntVec)}, shunts_by_bus, Int(sh.busIdx)), sh)
  end
  # 1) status contradiction: measured flow over an OPEN element
  n_checked_branches = 0
  for (k, br) in enumerate(net.branchVec)
    st = _branch_terminal_state(br)
    n_checked_branches += 1
    st == :closed && continue
    for m in get(flows_by_branch, k, Measurement[])
      ratio = abs(m.value) / max(m.sigma, eps())
      ratio > k_open || continue
      push!(findings, (stage = :precheck, kind = :open_element_with_flow, location = string("branch ", k, " (", busname(br.fromBus), "-", busname(br.toBus), ", ", st, ")"), evidence = string(m.id, " = ", round(m.value; sigdigits = 4), " at ", round(ratio; digits = 1), " sigma"), severity = _topology_severity(ratio, k_open)))
    end
  end
  n_checked_links = 0
  for (k, l) in enumerate(net.linkVec)
    n_checked_links += 1
    l.status == 1 && continue
    for m in get(flows_by_link, k, Measurement[])
      ratio = abs(m.value) / max(m.sigma, eps())
      ratio > k_open || continue
      push!(findings, (stage = :precheck, kind = :open_element_with_flow, location = string("link ", k, " (", busname(l.fromBus), "-", busname(l.toBus), ", open)"), evidence = string(m.id, " = ", round(m.value; sigdigits = 4), " at ", round(ratio; digits = 1), " sigma"), severity = _topology_severity(ratio, k_open)))
    end
  end

  # station -> does it carry any measured flow above the dead threshold?
  # Counted rather than flagged, split by the branch the row belongs to, so
  # check 2 can ask "live through an element OTHER than k" by subtraction
  # instead of rescanning every row per branch.
  station_live = Dict{Int,Bool}()
  live_by_station = Dict{Int,Int}()
  live_by_station_branch = Dict{Tuple{Int,Int},Int}()
  for m in active
    isflow(m) || continue
    abs(m.value) > k_dead * m.sigma || continue
    st = reps[endbus(m)]
    station_live[st] = true
    live_by_station[st] = get(live_by_station, st, 0) + 1
    key = (st, Int(m.branchIdx))
    live_by_station_branch[key] = get(live_by_station_branch, key, 0) + 1
  end
  # live flow rows at station `st` that do NOT belong to branch `k`
  live_elsewhere(st, k) = get(live_by_station, st, 0) - get(live_by_station_branch, (st, k), 0) > 0

  # 2) closed branch whose measurements all read dead while a terminal
  # station is otherwise live (needs at least one measurement per end;
  # an unmeasured end means no verdict, never a guess)
  for (k, br) in enumerate(net.branchVec)
    _branch_terminal_state(br) == :closed || continue
    rows = get(flows_by_branch, k, Measurement[])
    isempty(rows) && continue
    fb = Int(br.fromBus)
    tb = Int(br.toBus)
    has_from = any(endbus(m) == fb for m in rows)
    has_to = any(endbus(m) == tb for m in rows)
    (has_from && has_to) || continue
    all(abs(m.value) <= k_dead * m.sigma for m in rows) || continue
    # a terminal station must be live through OTHER elements, else the
    # whole neighbourhood is legitimately unloaded (false-positive guard)
    live_neighbour = live_elsewhere(reps[fb], k) || live_elsewhere(reps[tb], k)
    live_neighbour || continue
    push!(findings, (stage = :precheck, kind = :closed_element_without_flow, location = string("branch ", k, " (", busname(fb), "-", busname(tb), ", closed)"), evidence = string(length(rows), " flow/current row(s) all below ", k_dead, " sigma while a terminal station carries load"), severity = :warning))
  end

  # 3) closed link with voltage-magnitude measurements disagreeing across it
  vm_by_bus = Dict{Int,Measurement}()
  for m in active
    (m.typ == VmMeas && m.busIdx !== nothing) || continue
    haskey(vm_by_bus, m.busIdx) || (vm_by_bus[m.busIdx] = m)
  end
  for (k, l) in enumerate(net.linkVec)
    l.status == 1 || continue
    a = get(vm_by_bus, Int(l.fromBus), nothing)
    b = get(vm_by_bus, Int(l.toBus), nothing)
    (a === nothing || b === nothing) && continue
    σ = sqrt(a.sigma^2 + b.sigma^2)
    ratio = abs(a.value - b.value) / max(σ, eps())
    ratio > k_v || continue
    push!(findings, (stage = :precheck, kind = :closed_link_voltage_mismatch, location = string("link ", k, " (", busname(l.fromBus), "-", busname(l.toBus), ", closed)"), evidence = string("|", round(a.value; digits = 5), " - ", round(b.value; digits = 5), "| pu at ", round(ratio; digits = 1), " sigma"), severity = _topology_severity(ratio, k_v)))
  end

  # 4) node balance at COMPLETELY measured nodes (injection present, every
  # closed attached branch measured at this end); link-cluster members are
  # skipped (their balance lives on the fused station, stage 1 never
  # aggregates it)
  for (injT, flowT) in ((PinjMeas, PflowMeas), (QinjMeas, QflowMeas))
    for m in active
      (m.typ == injT && m.busIdx !== nothing) || continue
      b = m.busIdx
      length(get(members, reps[b], [b])) > 1 && continue
      complete = true
      flows = Measurement[]
      for k in get(closed_branches_by_bus, b, Int[])
        row = get(first_flow_row, (k, b, flowT), 0)
        if row == 0
          complete = false
          break
        end
        push!(flows, active[row])
      end
      (complete && !isempty(flows)) || continue
      # shunt contribution from the MEASURED voltage (the derived-shunt
      # principle): inj + y_sh * Vm^2 * base = sum(flows). A shunt bus
      # without a voltage measurement gets no verdict, never a guess.
      shTerm = 0.0
      shSkip = false
      for sh in get(shunts_by_bus, b, eltype(net.shuntVec)[])
        vmRow = get(vm_by_bus, b, nothing)
        if vmRow === nothing
          shSkip = true
          break
        end
        shTerm += (injT == PinjMeas ? real(sh.y_pu_shunt) : imag(sh.y_pu_shunt)) * vmRow.value^2 * net.baseMVA
      end
      shSkip && continue
      mismatch = abs(m.value + shTerm - sum(f.value for f in flows))
      σ = sqrt(m.sigma^2 + sum(f.sigma^2 for f in flows))
      ratio = mismatch / max(σ, eps())
      ratio > k_kcl || continue
      push!(findings, (stage = :precheck, kind = :kcl_violation, location = string("bus ", busname(b)), evidence = string(injT == PinjMeas ? "P" : "Q", " balance misses by ", round(mismatch; sigdigits = 4), " (", round(ratio; digits = 1), " sigma, ", length(flows), " measured branches)"), severity = _topology_severity(ratio, k_kcl)))
    end
  end

  summary = isempty(findings) ? string("topology precheck: no findings (", n_checked_branches, " branches, ", n_checked_links, " links checked)") : string("topology precheck: ", length(findings), " finding(s): ", join(("$(f.kind) at $(f.location)" for f in findings), "; "))
  return (findings = findings, summary = summary, n_checked_branches = n_checked_branches, n_checked_links = n_checked_links)
end

## candidate list for stage 3: explicit element references, or :auto from
## the stage-1 status-contradiction findings plus the switchable elements
## of stage-2 suspected stations (stage-2 stations first, open design
## decision 3). Elements are (kind = :branch | :link, idx).
function _topology_auto_candidates(net::Net, measurements::Vector{Measurement}; max_candidates::Int = 5, maxIte::Int = state_estimation_config().max_iter, tol::Float64 = state_estimation_config().tol)
  out = NamedTuple{(:kind, :idx),Tuple{Symbol,Int}}[]
  seen = Set{Tuple{Symbol,Int}}()
  addc(kind, idx) = begin
    (kind, idx) in seen && return nothing
    push!(seen, (kind, idx))
    push!(out, (kind = kind, idx = idx))
    return nothing
  end
  # stage-2 stations first: run the diagnostics fingerprint (with the
  # caller's solver settings, the tight defaults would under-iterate) and
  # collect every switchable element touching a suspected station
  diag = runse_diagnostics(net, measurements; maxIte = maxIte, tol = tol)
  reps, _ = _topology_station_map(net)
  for f in something(diag.topology_findings, NamedTuple[])
    f.kind == :topology_error_suspected_at_station || continue
    rep = f.station
    for (k, br) in enumerate(net.branchVec)
      (reps[Int(br.fromBus)] == rep || reps[Int(br.toBus)] == rep) && addc(:branch, k)
    end
    for (k, l) in enumerate(net.linkVec)
      (reps[Int(l.fromBus)] == rep || reps[Int(l.toBus)] == rep) && addc(:link, k)
    end
  end
  # then the stage-1 status contradictions (single elements)
  pre = validate_topology(net, measurements)
  for f in pre.findings
    f.kind in (:open_element_with_flow, :closed_element_without_flow) || continue
    mloc = match(r"^(branch|link) (\d+)", f.location)
    mloc === nothing && continue
    addc(Symbol(mloc.captures[1]), parse(Int, mloc.captures[2]))
  end
  length(out) > max_candidates && (out = out[1:max_candidates])
  return out
end

"""
    test_topology_hypotheses(net, measurements = net.measurements; candidates = :auto, max_candidates = 5, maxIte = ..., tol = ...) -> NamedTuple

Stage-3 topology hypothesis test: for each candidate element the
service-state is TOGGLED on a working copy of the net, the estimation is
re-run, and the chi-square band verdict before versus after decides
whether the hypothesis explains the data. The input net is never touched
(bitwise, test-enforced); the result is a RANKED RECOMMENDATION LIST and
nothing is ever switched automatically.

`candidates = :auto` collects the stage-2 suspected stations' switchable
elements first, then the stage-1 status contradictions, capped at
`max_candidates` (cost guard; per-candidate timing is reported). An
explicit vector of `(kind = :branch | :link, idx)` tuples is also
accepted. A hypothesis is `:hypothesis_supported` only when the toggled
run lands inside the Wilson-Hilferty band while the original run failed it
`:high`; several supported hypotheses are flagged `ambiguous`.
"""
function test_topology_hypotheses(net::Net, measurements::Vector{Measurement} = Measurement[m for m in net.measurements]; candidates = :auto, max_candidates::Int = 5, maxIte::Int = state_estimation_config().max_iter, tol::Float64 = state_estimation_config().tol)
  cands = candidates === :auto ? _topology_auto_candidates(net, measurements; max_candidates = max_candidates, maxIte = maxIte, tol = tol) : [(kind = Symbol(c.kind), idx = Int(c.idx)) for c in candidates]
  length(cands) > max_candidates && (cands = cands[1:max_candidates])
  # baseline on an untouched working copy (the caller's net stays pristine)
  base_net = deepcopy(net)
  base = runse!(base_net, measurements; maxIte = maxIte, tol = tol, updateNet = false, topologyPrecheck = false)
  base_verdict = _band_test_verdict(base.objectiveJ, base.dof)
  nby = _bus_name_by_idx(net)
  busname(b) = get(nby, Int(b), string(Int(b)))
  rows = NamedTuple[]
  for c in cands
    t0 = time_ns()
    work = deepcopy(net)
    label = ""
    current = :closed
    if c.kind == :branch
      (1 <= c.idx <= length(work.branchVec)) || continue
      br = work.branchVec[c.idx]
      current = _branch_terminal_state(br)
      setBranchStatus!(br, current != :closed)
      label = string("branch ", c.idx, " (", busname(net.branchVec[c.idx].fromBus), "-", busname(net.branchVec[c.idx].toBus), ")")
    else
      (1 <= c.idx <= length(work.linkVec)) || continue
      l = work.linkVec[c.idx]
      current = l.status == 1 ? :closed : :open
      l.status = current == :closed ? 0 : 1
      label = string("link ", c.idx, " (", busname(net.linkVec[c.idx].fromBus), "-", busname(net.linkVec[c.idx].toBus), ")")
    end
    hypothesis = current == :closed ? "actually open" : "actually closed"
    res = try
      runse!(work, measurements; maxIte = maxIte, tol = tol, updateNet = false, topologyPrecheck = false)
    catch err
      push!(rows, (element = label, kind = c.kind, idx = c.idx, current_status = current, hypothesis = hypothesis, j_before = base.objectiveJ, j_after = NaN, z_before = base_verdict.z_wh, z_after = NaN, verdict = :inconclusive, note = sprint(showerror, err), elapsed_s = (time_ns() - t0) / 1e9))
      continue
    end
    v = _band_test_verdict(res.objectiveJ, res.dof)
    # supported = the :high failure of the original run is GONE under the
    # hypothesis: inside the band, or below it (:low). A perfect hypothesis
    # on a noise-free set lands at J near zero, which the band flags :low;
    # for THIS question that is success, not doubt.
    verdict = if !res.converged
      :inconclusive
    elseif base_verdict.reason == :high && (v.passed || v.reason == :low)
      :hypothesis_supported
    else
      :not_supported
    end
    push!(rows, (element = label, kind = c.kind, idx = c.idx, current_status = current, hypothesis = hypothesis, j_before = base.objectiveJ, j_after = res.objectiveJ, z_before = base_verdict.z_wh, z_after = v.z_wh, verdict = verdict, note = "", elapsed_s = (time_ns() - t0) / 1e9))
  end
  sort!(rows; by = r -> (r.verdict == :hypothesis_supported ? 0 : 1, isnan(r.j_after) ? Inf : r.j_after))
  supported = count(r -> r.verdict == :hypothesis_supported, rows)
  return (recommendations = rows, base_j = base.objectiveJ, base_dof = base.dof, base_verdict = base_verdict.reason, ambiguous = supported > 1, n_candidates = length(rows))
end
