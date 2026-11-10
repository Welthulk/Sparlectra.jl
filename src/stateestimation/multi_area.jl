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

# file: src/stateestimation/multi_area.jl
# purpose: multi-area state estimation, two-level method (1.2.0): the
#          partition of a net into areas, the extended subnet of each area,
#          the area run with the foreign ends of the tie branches held fixed,
#          the border information (Schur complement of the area's gain
#          matrix onto the border columns) and the coordinator that solves
#          the border states from the sum of the area models.

# ---------------------------------------------------------------------------
# Partition
# ---------------------------------------------------------------------------

"""
    AreaPartition

The areas of a net for the multi-area state estimation.

- `area_of_bus`: area per bus (bus numbering of the net).
- `areas`: the area numbers, ascending.
- `area_names`: the name per entry of `areas` (`net.area_names`: the
  operator model or region of a CGMES import, SCF `sparlectra.areas`);
  empty strings without names and with an area list.
- `tie_branches`: closed branches whose ends lie in two areas.
- `border_buses`: the ends of the tie branches, ascending.
- `notes`: what the partition changed, e.g. a link group the area list cut,
  moved to the area of its representative.
"""
struct AreaPartition
  area_of_bus::Vector{Int}
  areas::Vector{Int}
  area_names::Vector{String}
  tie_branches::Vector{Int}
  border_buses::Vector{Int}
  notes::Vector{String}
end

"""
    area_partition(net; areas = nothing) -> AreaPartition

Areas of `net` for the multi-area state estimation. Without `areas` every
bus takes the area stored with it (`area` of the MATPOWER bus table, SCF
`extra.area`, the `area` keyword of `addBus!`); a bus without one is an
error that names it. `areas` is a map bus name => area number and must name
every bus.

A group of buses joined by closed links is never cut: when the list puts its
members in different areas, the group goes to the area of its representative
(the lowest bus index), with a note. Open links and one-sided open branches
do not couple areas.
"""
function area_partition(net::Net; areas::Union{Nothing,AbstractDict} = nothing)::AreaPartition
  n = length(net.nodeVec)
  area = zeros(Int, n)
  missingNames = String[]
  for (b, nd) in enumerate(net.nodeVec)
    name = getCompName(nd.comp)
    a = areas === nothing ? nd._area : get(areas, name, nothing)
    if a === nothing
      push!(missingNames, name)
    else
      area[b] = Int(a)
    end
  end
  if !isempty(missingNames)
    src = areas === nothing ? "no stored area" : "not in the area list"
    shown = join(first(missingNames, 5), ", ")
    throw(ArgumentError("area_partition: $(length(missingNames)) bus(es) with $(src): $(shown)$(length(missingNames) > 5 ? ", ..." : "")"))
  end
  notes = String[]
  reps = _active_link_representative_map(net)
  for b = 1:n
    rb = reps[b]
    if area[b] != area[rb]
      push!(notes, "bus $(getCompName(net.nodeVec[b].comp)) moved from area $(area[b]) to area $(area[rb]): a closed link joins it to $(getCompName(net.nodeVec[rb].comp))")
      area[b] = area[rb]
    end
  end
  ties = Int[]
  for (k, br) in enumerate(net.branchVec)
    _branch_terminal_state(br) === :closed || continue
    area[Int(br.fromBus)] != area[Int(br.toBus)] && push!(ties, k)
  end
  border = sort!(unique!(vcat(Int[Int(net.branchVec[k].fromBus) for k in ties], Int[Int(net.branchVec[k].toBus) for k in ties])))
  nums = sort!(unique(area))
  # the names belong to the numbers the net stores; a list numbers its own
  names = String[areas === nothing ? get(net.area_names, a, "") : "" for a in nums]
  return AreaPartition(area, nums, names, ties, border, notes)
end

## The area that owns a measurement row: a bus row
## the area of its bus, a branch row the area of its measured end, a link
## row the area of its link group (never cut, see `area_partition`).
function _se_area_owner(net::Net, part::AreaPartition, m::Measurement)::Int
  if m.linkIdx !== nothing
    return part.area_of_bus[Int(net.linkVec[m.linkIdx].fromBus)]
  elseif m.branchIdx !== nothing
    br = net.branchVec[m.branchIdx]
    return part.area_of_bus[m.direction === :to ? Int(br.toBus) : Int(br.fromBus)]
  elseif m.busIdx !== nothing
    return part.area_of_bus[m.busIdx]
  end
  return 0
end

# ---------------------------------------------------------------------------
# Extended subnet of one area
# ---------------------------------------------------------------------------

## The extended subnet of area `a`: its own buses, then the copies of the
## foreign ends of its tie branches (decision E3). A copy carries no
## prosumer, shunt or link: only the tie branch reaches it, and no row of
## the area is a bus row of the copy. Branches: the area's own (both ends,
## or the closed end of a one-sided open branch, in the area) and its tie
## branches. `V` gives the start voltage per bus of the net; every border
## bus of the subnet (own border buses and copies) keeps its value for the
## whole run (`pinned`): the area estimates its interior for given border
## values, the reduced problem of the two-level method. Without the slack
## of the net among its own buses, the first copy becomes the angle
## reference of the subnet, at its value from `V`.
function _se_area_subnet(net::Net, part::AreaPartition, a::Int, V::Vector{ComplexF64}; base::Union{Nothing,Net} = nothing)
  own = [b for b in eachindex(part.area_of_bus) if part.area_of_bus[b] == a]
  ownset = Set(own)
  ties = [k for k in part.tie_branches if Int(net.branchVec[k].fromBus) in ownset || Int(net.branchVec[k].toBus) in ownset]
  copies = sort!(unique!([Int(net.branchVec[k].fromBus) in ownset ? Int(net.branchVec[k].toBus) : Int(net.branchVec[k].fromBus) for k in ties]))
  inner = Int[]
  for (k, br) in enumerate(net.branchVec)
    st = _branch_terminal_state(br)
    f, t = Int(br.fromBus), Int(br.toBus)
    keep = st === :closed ? (f in ownset && t in ownset) : st === :open_to ? f in ownset : st === :open_from ? t in ownset : false
    keep && push!(inner, k)
  end
  buses = vcat(own, copies)
  branches = sort!(vcat(inner, ties))
  slackOwn = any(b -> getNodeType(net.nodeVec[b]) == Slack, own)
  row = (buses = buses, branches = branches, n_bus = length(buses), n_ref = slackOwn ? 1 : 0, chosen_ref_bus = 0)
  inet = _prepare_island_net(net, row; base = base)
  nOwn = length(own)
  copyIdx = collect((nOwn+1):length(buses))
  copyset = Set(copyIdx)
  # a copy is a plain bus: no unit (a foreign slack would become a second
  # reference), no shunt
  inet.prosumpsVec = [ps for ps in inet.prosumpsVec if !(Int(ps.comp.cFrom_bus) in copyset)]
  inet.shuntVec = [sh for sh in inet.shuntVec if !(Int(sh.busIdx) in copyset)]
  empty!(inet.shuntDict)
  for (k, sh) in enumerate(inet.shuntVec)
    inet.shuntDict[Int(sh.busIdx)] = k
  end
  inet.slackVec = [s for s in inet.slackVec if !(s in copyset)]
  busmap = Dict{Int,Int}(old => new for (new, old) in enumerate(buses))
  # links: only those inside the area (a closed group is never cut; an open
  # link across the border does not couple and is left out)
  linkmap = Dict{Int,Int}()
  links = eltype(net.linkVec)[]
  for (k, l) in enumerate(net.linkVec)
    (Int(l.fromBus) in ownset && Int(l.toBus) in ownset) || continue
    nl = deepcopy(l)
    nl.fromBus = busmap[Int(l.fromBus)]
    nl.toBus = busmap[Int(l.toBus)]
    nl.linkIdx = length(links) + 1
    push!(links, nl)
    linkmap[k] = nl.linkIdx
  end
  inet.linkVec = links
  if !slackOwn
    isempty(copyIdx) && throw(ArgumentError("multi-area estimation: area $(a) has neither the slack nor a tie branch"))
    ref = first(copyIdx)
    refname = getCompName(inet.nodeVec[ref].comp)
    addProsumer!(net = inet, busName = refname, type = "EXTERNALNETWORKINJECTION", referencePri = refname, vm_pu = abs(V[buses[ref]]), va_deg = rad2deg(angle(V[buses[ref]])))
  end
  refreshBusTypesFromProsumers!(inet)
  slackOwn || push!(inet.slackVec, first(copyIdx))
  for (k, b) in enumerate(buses)
    inet.nodeVec[k]._vm_pu = abs(V[b])
    inet.nodeVec[k]._va_deg = rad2deg(angle(V[b]))
  end
  branchmap = Dict{Int,Int}(old => new for (new, old) in enumerate(branches))
  borderset = Set(part.border_buses)
  pinned = vcat([k for (k, b) in enumerate(own) if b in borderset], copyIdx)
  return (inet = inet, buses = buses, own = own, copies = copies, pinned = pinned, busmap = busmap, branchmap = branchmap, linkmap = linkmap, shuntpos = Int[])
end

# ---------------------------------------------------------------------------
# Area run and border information
# ---------------------------------------------------------------------------

## The border model of one area from the last linearization of its run:
## per border variable (bus of the net, :va or :vm) the value at the area's
## point, the reduced gradient and the Schur complement of the gain matrix
## onto the border columns (the interior eliminated). Angles in rad,
## absolute (slack frame of the net), magnitudes in pu. With the gradient g
## of J/2 and the gain G = H'WH: S = G_BB - G_BI G_II^-1 G_IB,
## g_red = g_B - G_BI G_II^-1 g_I.
##
## The angle column of the subnet's reference is not a state of the run.
## Moving that bus by d with every other absolute angle held moves the
## reference by d and every relative angle by -d; the angle rows (Va, Ia)
## read absolute angles and stay, so the column is the sum over those rows'
## unit entries minus the sum of all other angle columns.
function _se_area_border_model(lin, gbus::Vector{Int}, isborder::Function)
  cap = lin.cap
  cap === nothing && error("multi-area estimation: the area run returned no linearization")
  H = cap.H
  w = cap.w
  r = cap.r
  x = lin.x
  nb = lin.nbus
  s = lin.slackIdx
  nθ = nb - 1
  refRad = deg2rad(lin.vaRefDeg)
  vacol(i) = i < s ? i : i - 1
  ncol = size(H, 2)
  # absolute-angle rows: their entry of the reference column
  absRow = Float64[(m.typ == VaMeas || m.typ == IaMeas) ? 1.0 : 0.0 for m in lin.meas]
  borderVars = Tuple{Int,Symbol}[]
  borderVals = Float64[]
  bcols = Int[]          # state columns of the border variables (0: the reference angle)
  isBorderCol = falses(ncol)
  for i = 1:nb
    gb = gbus[lin.buses[i]]
    isborder(gb) || continue
    push!(borderVars, (gb, :va))
    push!(borderVals, i == s ? refRad : x[vacol(i)] + refRad)
    if i == s
      push!(bcols, 0)
    else
      push!(bcols, vacol(i))
      isBorderCol[vacol(i)] = true
    end
    push!(borderVars, (gb, :vm))
    push!(borderVals, x[nθ+i])
    push!(bcols, nθ + i)
    isBorderCol[nθ+i] = true
  end
  icols = [c for c = 1:ncol if !isBorderCol[c]]
  HBd = Matrix{Float64}(undef, size(H, 1), length(bcols))
  for (k, c) in enumerate(bcols)
    if c == 0
      col = copy(absRow)
      for i = 1:nb
        i == s && continue
        col .-= Vector(H[:, vacol(i)])
      end
      HBd[:, k] = col
    else
      HBd[:, k] = Vector(H[:, c])
    end
  end
  HI = H[:, icols]
  Wd = Diagonal(w)
  GII = Matrix(HI' * (Wd * HI))
  GIB = Matrix(HI' * (Wd * HBd))
  GBB = HBd' * (Wd * HBd)
  gI = -(HI' * (w .* r))
  gB = -(HBd' * (w .* r))
  F = nothing
  if isempty(icols)
    S = GBB
    g = gB
  else
    F = cholesky(Symmetric(GII); check = false)
    issuccess(F) || (F = lu(GII))
    X = F \ GIB
    y = F \ gI
    S = GBB - GIB' * X
    g = gB - GIB' * y
  end
  # the blocks stay for the residual covariance of the area's rows
  # (`_se_multi_area_residuals`)
  return (vars = borderVars, vals = borderVals, S = Symmetric(0.5 .* (S .+ S')), g = g, ninterior = length(icols), F = F, GIB = GIB, HI = HI, HB = HBd, w = w, r = r)
end

## The area runs color the finite-difference Jacobian from this state count
## on (the estimator itself from 1000): an area of sp_case1354 has about 560
## states and 2000 rows, column by column its Jacobian took most of the run
## (2.7 s of the passes against 0.1 s for runse! on the whole net, which is
## colored). The colored Jacobian is the same bit for bit.
const _SE_AREA_FD_COLORING_MIN_STATES = 200

## The tap exclusions of a multi-area run (concept 3.3, rule 1): those of
## the configuration plus every released transformer on a tie branch
## (reason `:tie_branch`; its step would couple two areas). `ties`: the tie
## branches as indices of `net` (0 = not in this net).
function _se_area_tap_exclusions(net::Net, cfg::StateEstimationConfig, ties::Vector{Int})::Dict{String,Symbol}
  excl = _tap_exclusion_map(net, cfg.tap_estimation)
  cfg.tap_estimation.enabled || return excl
  for k in ties
    k == 0 && continue
    br = net.branchVec[k]
    (br.ratio != 0.0 && br.tap_est_mode !== :none) && (excl[br.comp.cID] = :tie_branch)
  end
  return excl
end

## One area run: the extended subnet with the copies at `V`, the area's own
## rows, the copies pinned, start from `V` (flatstart off). Returns the
## voltages of the own buses (net numbering) and the border model.
function _se_area_run(net::Net, part::AreaPartition, a::Int, rows::Vector{Int}, measurements::Vector{Measurement}, cfg::StateEstimationConfig, V::Vector{ComplexF64}; base::Union{Nothing,Net} = nothing, tap_fix::Bool = false)
  sub = _se_area_subnet(net, part, a, V; base = base)
  meas, kept = _se_remap_island_measurements(sub, measurements, rows)
  # update_taps on the fixation run only: it writes the fixed steps into the
  # subnet, from where `_se_area_fix_taps!` takes them
  acfg = _copy_se_config_with(cfg; flatstart = false, update_net = true, update_taps = tap_fix)
  lin = Ref{Any}(nothing)
  excl = _se_area_tap_exclusions(sub.inet, cfg, [get(sub.branchmap, k, 0) for k in part.tie_branches])
  res = _runse_with_config!(sub.inet, meas, acfg; machine_rows = false, tap_excl = excl, pinned_buses = sub.pinned, linearization = lin, tap_fix = tap_fix, fd_coloring_min_states = _SE_AREA_FD_COLORING_MIN_STATES)
  Vown = Dict{Int,ComplexF64}()
  for (k, b) in enumerate(sub.own)
    nd = sub.inet.nodeVec[k]
    Vown[b] = nd._vm_pu * cis(deg2rad(nd._va_deg))
  end
  borderset = Set(part.border_buses)
  model = _se_area_border_model(lin[], sub.buses, b -> b in borderset)
  return (V = Vown, model = model, result = res, sub = sub, lin = lin[], kept = kept)
end

# ---------------------------------------------------------------------------
# Coordinator
# ---------------------------------------------------------------------------

## The area runs of one pass (decision E7): areas of at least
## `power_flow.islands.parallel_min_buses` buses run on their own tasks when
## two or more qualify, Julia has more than one thread and
## `runtime.parallel.enabled` is on, as for the islands of `runse!`; the rest
## follow on this task. Every run cuts its own subnet from the shared base and
## reads `V` and the rows only; nothing mutable is shared. The serial path is
## the same function.
function _se_area_runs(net::Net, part::AreaPartition, rowsOf::Dict{Int,Vector{Int}}, measurements::Vector{Measurement}, cfg::StateEstimationConfig, V::Vector{ComplexF64}, base::Net; tap_fix::Bool = false)::Dict{Int,Any}
  runs = Dict{Int,Any}()
  active = active_sparlectra_config()
  pr = active.runtime.parallel
  min_buses = active.powerflow.islands_parallel_min_buses
  nbuses = Dict{Int,Int}(a => count(==(a), part.area_of_bus) for a in part.areas)
  big = [a for a in part.areas if nbuses[a] >= min_buses]
  cap = parallel_max_tasks(pr)
  if pr.enabled && Threads.nthreads() > 1 && cap > 1 && length(big) >= 2
    tasks = Dict{Int,Task}()
    sem = Base.Semaphore(max(1, cap))
    for a in sort(big; by = x -> nbuses[x], rev = true)
      let area = a
        tasks[area] = Threads.@spawn begin
          Base.acquire(sem)
          try
            _se_area_run(net, part, area, rowsOf[area], measurements, cfg, V; base = base, tap_fix = tap_fix)
          finally
            Base.release(sem)
          end
        end
      end
    end
    for a in big
      runs[a] = fetch(tasks[a])
    end
  end
  for a in part.areas
    haskey(runs, a) || (runs[a] = _se_area_run(net, part, a, rowsOf[a], measurements, cfg, V; base = base, tap_fix = tap_fix))
  end
  return runs
end

## The fixation of the released regulators after the continuous passes:
## each area with an estimated regulator runs once more at the converged
## voltages `V` with the fixation of runse! (`tap_fix = true`,
## update_taps), and the fixed steps are copied from its subnet into
## `work` (the branch fields, never the bus indices), where the regulator is
## no longer released. Returns the tap rows of these runs with their area.
function _se_area_fix_taps!(work::Net, net::Net, part::AreaPartition, rowsOf::Dict{Int,Vector{Int}}, measurements::Vector{Measurement}, cfg::StateEstimationConfig, V::Vector{ComplexF64}, base::Net)
  rows = NamedTuple[]
  seen = Set{Tuple{Int,Int}}()
  keepFields = (:branchIdx, :fromBus, :toBus, :comp)
  for a in part.areas
    run = _se_area_run(net, part, a, rowsOf[a], measurements, cfg, V; base = base, tap_fix = true)
    te = run.result.tapEstimates
    te === nothing && continue
    back = Dict{Int,Int}(new => old for (old, new) in run.sub.branchmap)
    for t in te
      # a tie transformer sits in both its areas' subnets: one row
      key = (back[t.branch], t.regulator)
      key in seen || (push!(seen, key); push!(rows, merge(t, (area = a, branch = back[t.branch]))))
      (t.status === :estimated && t.fixed) || continue
      k = back[t.branch]
      src = run.sub.inet.branchVec[t.branch]
      dst = work.branchVec[k]
      for fld in fieldnames(typeof(dst))
        fld in keepFields && continue
        setfield!(dst, fld, deepcopy(getfield(src, fld)))
      end
      setTapEstimation!(work; trafo = k, enabled = false)
    end
  end
  return rows
end

## The tap rows of area runs in the numbering of the net, with their area.
## A tie transformer appears in the subnets of both its areas (reason
## `:tie_branch` in each); its row is kept once, from the first area.
function _se_area_tap_rows(runs::Dict{Int,Any}, part::AreaPartition)::Vector{NamedTuple}
  rows = NamedTuple[]
  seen = Set{Tuple{Int,Int}}()
  for a in part.areas
    te = runs[a].result.tapEstimates
    te === nothing && continue
    back = Dict{Int,Int}(new => old for (old, new) in runs[a].sub.branchmap)
    for t in te
      key = (back[t.branch], t.regulator)
      key in seen && continue
      push!(seen, key)
      push!(rows, merge(t, (area = a, branch = back[t.branch])))
    end
  end
  return rows
end

## The release of the tap regulators per area against the decision of the
## integrated estimate on the same rows (each
## area decides, the run compares and names a difference). One row per
## regulator of either side: `branch`, `regulator`, `name`, `changer`,
## `area`, `multi_area` and `integrated` as (status, reason), `same`.
function _se_tap_release_compare(maRows, intRows)::Vector{NamedTuple}
  key(t) = (t.branch, t.regulator)
  ma = Dict{Tuple{Int,Int},NamedTuple}()
  for t in something(maRows, NamedTuple[])
    haskey(ma, key(t)) || (ma[key(t)] = t)
  end
  it = Dict{Tuple{Int,Int},NamedTuple}(key(t) => t for t in something(intRows, NamedTuple[]))
  out = NamedTuple[]
  for k in sort!(collect(union(keys(ma), keys(it))))
    a = get(ma, k, nothing)
    b = get(it, k, nothing)
    t = something(a, b)
    sa = a === nothing ? (:absent, :none) : (a.status, a.reason)
    sb = b === nothing ? (:absent, :none) : (b.status, b.reason)
    push!(out, (branch = k[1], regulator = k[2], name = t.name, changer = t.changer, area = a === nothing ? 0 : a.area, multi_area = sa, integrated = sb, same = sa == sb))
  end
  return out
end

## The coordinator's linear system from the area models: sum_a S_a b =
## sum_a (S_a b_a - g_a) over the border variables `vix` (variable =>
## position); a model entry outside `vix` (the slack angle) is dropped.
function _se_coordinator_system(models::Vector, vix::Dict{Tuple{Int,Symbol},Int}, nv::Int)
  A = zeros(nv, nv)
  rhs = zeros(nv)
  for md in models
    idx = [get(vix, v, 0) for v in md.vars]
    for (p, ip) in enumerate(idx)
      ip == 0 && continue
      rhs[ip] -= md.g[p]
      for (q, iq) in enumerate(idx)
        iq == 0 && continue
        A[ip, iq] += md.S[p, q]
        rhs[ip] += md.S[p, q] * md.vals[q]
      end
    end
  end
  return A, rhs
end

## The two-level estimation (concept 1.2.0, section 3.2). Per pass every
## area estimates its interior with its own rows and all border buses of its
## subnet held at the coordinator's values b; the coordinator then takes one
## Gauss-Newton step on all border variables from the sum of the area
## models, sum_a S_a d = -sum_a g_a: Newton on the reduced objective
## sum_a min_interior J_a. At the fixed point the reduced gradients sum to
## zero, the optimality condition of the integrated estimate. Stops when no
## border variable moves by more than `tol`. (Holding only the copies and
## leaving the own border free gave a model at the area's own optimum, not
## at b; with noise the passes stalled at 2e-4.)
## `check`: refuse an area whose rows do not determine its interior.
function _runse_two_level(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig, part::AreaPartition; max_passes::Int = 20, V0::Union{Nothing,Vector{ComplexF64}} = nothing, check::Bool = true)
  nbus = length(net.nodeVec)
  slack = findall(nd -> getNodeType(nd) == Slack, net.nodeVec)
  length(slack) == 1 || throw(ArgumentError("multi-area estimation: the net needs exactly one slack bus, found $(length(slack))"))
  slackBus = only(slack)
  cfg.zero_injection === :constraint && throw(ArgumentError("multi-area estimation: state_estimation.zero_injection = constraint is not supported, use on or off"))
  rowsOf = Dict{Int,Vector{Int}}(a => Int[] for a in part.areas)
  for (k, m) in enumerate(measurements)
    m.active || continue
    a = _se_area_owner(net, part, m)
    a == 0 && continue
    push!(rowsOf[a], k)
  end
  pmuAreas = [a for a in part.areas if any(k -> measurements[k].typ in (VaMeas, IaMeas), rowsOf[a])]
  length(pmuAreas) > 1 && throw(ArgumentError("multi-area estimation: PMU angle rows in more than one area ($(join(pmuAreas, ", "))) need a shared offset, not supported yet"))
  # start: the stored voltages, or flat at the slack angle
  vaRef = deg2rad(_se_stored_value(net.nodeVec[slackBus]._va_deg, 0.0))
  V = if V0 !== nothing
    copy(V0)
  elseif cfg.flatstart
    [ComplexF64(b == slackBus ? _se_stored_value(net.nodeVec[b]._vm_pu, 1.0) : 1.0, 0.0) * cis(vaRef) for b = 1:nbus]
  else
    [_se_stored_value(nd._vm_pu, 1.0) * cis(deg2rad(_se_stored_value(nd._va_deg, rad2deg(vaRef)))) for nd in net.nodeVec]
  end
  # border variables: angle (not of the slack) and magnitude per border bus
  vars = Tuple{Int,Symbol}[]
  for b in part.border_buses
    b == slackBus || push!(vars, (b, :va))
    push!(vars, (b, :vm))
  end
  vix = Dict{Tuple{Int,Symbol},Int}(v => k for (k, v) in enumerate(vars))
  for a in part.areas
    isempty(rowsOf[a]) && throw(ArgumentError("multi-area estimation: area $(a) has no active measurement"))
    if check
      msg = _se_area_check(net, part, a, rowsOf[a], measurements, V)
      msg === nothing || throw(ArgumentError("multi-area estimation: " * msg))
    end
  end
  # one copy of the net without the containers every subnet replaces (the
  # net does not change during the passes)
  base = _se_island_base(net)
  # phase 1: released regulators stay continuous interior states of their
  # area, so the area models carry them and the fixed point is the
  # continuous estimate of the integrated run. Fixing them per area inside
  # every pass drifted: on sp_case118 with nine regulators the passes
  # settled 1 to 5 steps and 0.03 pu away from runse!, since a
  # fixed step and the border magnitudes made up for each other.
  ph = _se_two_level_passes!(net, base, part, rowsOf, measurements, cfg, V, vars, vix, max_passes; tap_fix = false)
  history = copy(ph.history)
  passes = ph.passes
  converged = ph.converged
  runs = ph.runs
  work = net
  tapRows = NamedTuple[]
  if converged && any(a -> runs[a].result.tapEstimates !== nothing && any(t -> t.status === :estimated, runs[a].result.tapEstimates), part.areas)
    # phase 2: each area fixes its regulators once at the converged border
    # (the rule of runse!), the steps go into a copy of the net, and the
    # passes go on with them fixed, like the fixation run of runse!
    work = deepcopy(net)
    tapRows = _se_area_fix_taps!(work, net, part, rowsOf, measurements, cfg, V, base)
    ph2 = _se_two_level_passes!(work, _se_island_base(work), part, rowsOf, measurements, cfg, V, vars, vix, max(max_passes - passes, 1); tap_fix = true)
    append!(history, ph2.history)
    passes += ph2.passes
    converged = ph2.converged
    runs = ph2.runs
  elseif converged
    # no regulator estimated in any area: the rows of the continuous runs
    # still carry the release decision of each area
    tapRows = _se_area_tap_rows(runs, part)
  end
  return (V = V, converged = converged, passes = passes, history = history, runs = runs, rows = rowsOf, vars = vars, vix = vix, net = work, tap_rows = tapRows)
end

## The passes of the two-level method on `net` from the voltages `V`
## (updated in place), see `_runse_two_level`.
function _se_two_level_passes!(net::Net, base::Net, part::AreaPartition, rowsOf::Dict{Int,Vector{Int}}, measurements::Vector{Measurement}, cfg::StateEstimationConfig, V::Vector{ComplexF64}, vars::Vector{Tuple{Int,Symbol}}, vix::Dict{Tuple{Int,Symbol},Int}, max_passes::Int; tap_fix::Bool)
  nv = length(vars)
  passes = 0
  change = Inf
  history = Float64[]
  runs = nothing
  converged = false
  tol = max(cfg.tol, cfg.jac_eps)
  for pass = 1:max_passes
    passes = pass
    runs = _se_area_runs(net, part, rowsOf, measurements, cfg, V, base; tap_fix = tap_fix)
    A, rhs = _se_coordinator_system([runs[a].model for a in part.areas], vix, nv)
    # a border variable no area sees has an empty row: name the bus
    for (k, (b, q)) in enumerate(vars)
      A[k, k] > 0.0 || throw(ArgumentError("multi-area estimation: border bus $(getCompName(net.nodeVec[b].comp)) ($(q)) is determined by no area"))
    end
    F = cholesky(Symmetric(A); check = false)
    issuccess(F) || throw(ArgumentError("multi-area estimation: the border states are not determined by the areas together (coordinator matrix not positive definite)"))
    bnew = F \ rhs
    # interior and own border voltages from the areas, then the border
    # buses from the coordinator
    for a in part.areas
      for (b, v) in runs[a].V
        V[b] = v
      end
    end
    change = 0.0
    for (k, (b, q)) in enumerate(vars)
      vm, va = abs(V[b]), angle(V[b])
      if q === :va
        change = max(change, abs(bnew[k] - va))
        va = bnew[k]
      else
        change = max(change, abs(bnew[k] - vm))
        vm = bnew[k]
      end
      V[b] = vm * cis(va)
    end
    push!(history, change)
    if change < tol
      converged = true
      break
    end
  end
  return (converged = converged, passes = passes, history = history, runs = runs)
end

# ---------------------------------------------------------------------------
# Observability of an area for fixed border values (decision E5)
# ---------------------------------------------------------------------------

## Whether the own rows of area `a` determine its interior when the border
## buses of its subnet are given: the decoupled check of
## `evaluate_global_observability` on the subnet, with the border buses
## measured by pseudo rows (Vm and Va). An area that fails is refused with
## its deficits.
function _se_area_check(net::Net, part::AreaPartition, a::Int, rows::Vector{Int}, measurements::Vector{Measurement}, V::Vector{ComplexF64})
  sub = _se_area_subnet(net, part, a, V)
  meas, _ = _se_remap_island_measurements(sub, measurements, rows)
  for k in sub.pinned
    push!(meas, Measurement(typ = VmMeas, value = abs(V[sub.buses[k]]), sigma = 1e-3, busIdx = k, id = "border_vm_$(k)"))
    push!(meas, Measurement(typ = VaMeas, value = rad2deg(angle(V[sub.buses[k]])), sigma = 1e-3, busIdx = k, id = "border_va_$(k)"))
  end
  obs = evaluate_global_observability(sub.inet, meas; observable_islands = false)
  obs.numerical_observable && return nothing
  return "area $(a) is not observable for fixed border values (angle rank deficit $(obs.angle_rank_deficit), magnitude rank deficit $(obs.magnitude_rank_deficit), $(length(rows)) own row(s))"
end

# ---------------------------------------------------------------------------
# Result and entry point
# ---------------------------------------------------------------------------

"""
    MultiAreaSEResult

Result of [`runse_multi_area!`](@ref).

- `voltages`: estimated complex voltages per bus (bus numbering of the net,
  slack angle frame of the net).
- `converged`, `passes`: whether the border states settled (no change above
  `state_estimation.tol`) and after how many passes; `border_change` holds
  the largest change per pass.
- `objectiveJ`, `dof`: J over all rows (the areas' rows are disjoint) and the
  degrees of freedom of the whole set.
- `areas`: one row per area: `area`, `name` (see [`AreaPartition`](@ref)),
  `buses`, `border_buses`, `rows` (count
  of own active rows), `j` and `iterations` of its last run, `converged`,
  and `result`, the [`SEResult`](@ref) of that run on the extended subnet.
- `borders`: one row per border bus: `bus`, `name`, `area`, `vm_pu`,
  `va_deg`.
- `partition`: the [`AreaPartition`](@ref) used.
- `tapEstimates`: the released tap regulators, one row per changer as in
  `SEResult.tapEstimates` plus its `area` (`branch` in the numbering of the
  net); `nothing` without released regulators. The passes keep the
  regulators as continuous states; at the converged border each area fixes
  them with the rule of `runse!`, and the passes go on with the fixed steps.
- `comparison`: with `compare_integrated = true` the integrated estimate on
  the same rows: `rows` per state (`bus`, `name`, `quantity` (:vm in pu or
  :va in degrees), `multi_area`, `integrated`, `difference`, `sigma`, the
  standard deviation of the state in the integrated estimate, and
  `difference_sigma`), `max_dvm_pu`, `max_dva_deg`, `integrated_converged`,
  and the check of the tap release: `tap_release`, one row per released
  regulator (`branch`, `regulator`, `name`, `changer`, `area`, `multi_area`
  and `integrated` as (status, reason), `same`) and
  `tap_release_differences`, the number of rows that differ (named in a
  warning); `nothing` otherwise. Each area decides the release of its own
  regulators; without the comparison that decision is not checked.
"""
struct MultiAreaSEResult
  voltages::Vector{ComplexF64}
  converged::Bool
  passes::Int
  border_change::Vector{Float64}
  objectiveJ::Float64
  dof::Int
  areas::Vector{NamedTuple}
  borders::Vector{NamedTuple}
  partition::AreaPartition
  tapEstimates::Union{Nothing,Vector{NamedTuple}}
  comparison::Union{Nothing,NamedTuple}
end

## The partition of a multi-area run: the map or list file the caller
## passed, else `state_estimation.multi_area.partition` (`bus_area` or the
## path of an area list).
function _se_configured_partition(net::Net, cfg::StateEstimationConfig, areas)::AreaPartition
  areas isa AbstractDict && return area_partition(net; areas = areas)
  src = areas isa AbstractString ? String(areas) : cfg.multi_area.partition
  src == "bus_area" && return area_partition(net)
  return area_partition(net; areas = read_area_list(src))
end

"""
    read_area_list(path) -> Dict{String,Int}

Read an area list for [`area_partition`](@ref): one line per bus with the bus
name and the area number, separated by `;`, `,` or a tab. Empty lines, lines
starting with `#` and a header line (first data line whose second field is
not a number) are skipped. A malformed line or a bus named twice is an
`ArgumentError` with the line number.
"""
function read_area_list(path::AbstractString)::Dict{String,Int}
  isfile(path) || throw(ArgumentError("read_area_list: no such file: $(path)"))
  out = Dict{String,Int}()
  header = true
  for (ln, raw) in enumerate(eachline(path))
    line = strip(raw)
    (isempty(line) || startswith(line, "#")) && continue
    f = strip.(split(line, r"[;,\t]"))
    length(f) >= 2 || throw(ArgumentError("read_area_list: line $(ln) of $(path) needs a bus name and an area number: $(line)"))
    a = tryparse(Int, f[2])
    if a === nothing
      header && (header = false; continue)
      throw(ArgumentError("read_area_list: line $(ln) of $(path): area $(repr(f[2])) is not an integer"))
    end
    header = false
    haskey(out, f[1]) && throw(ArgumentError("read_area_list: line $(ln) of $(path): bus $(f[1]) is named twice"))
    out[String(f[1])] = a
  end
  return out
end

## The integrated estimate on the same rows for the comparison: one run of
## the estimator core on a copy of the net, with its last linearization for
## the standard deviations of the states.
function _se_multi_area_compare(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig, part::AreaPartition, V::Vector{ComplexF64}; tap_rows = NamedTuple[])
  ni = deepcopy(net)
  lin = Ref{Any}(nothing)
  res = _runse_with_config!(ni, measurements, _copy_se_config_with(cfg; update_net = true); machine_rows = false, tap_excl = _se_area_tap_exclusions(ni, cfg, part.tie_branches), linearization = lin)
  Vi = ComplexF64[_se_stored_value(nd._vm_pu, NaN) * cis(deg2rad(_se_stored_value(nd._va_deg, NaN))) for nd in ni.nodeVec]
  l = lin[]
  nb, s = l.nbus, l.slackIdx
  vacol(i) = i < s ? i : i - 1
  rowsUsed = [k for k in eachindex(l.cap.w) if l.cap.w[k] > 0.0]
  cols = Int[]
  for i = 1:nb
    i == s || push!(cols, vacol(i))
    push!(cols, nb - 1 + i)
  end
  sig = _se_state_sigmas(l.cap.H, l.cap.w, Int[], rowsUsed, cols)
  out = NamedTuple[]
  maxVm = 0.0
  maxVa = 0.0
  for i = 1:nb
    b = l.buses[i]
    name = getCompName(net.nodeVec[b].comp)
    dvm = abs(V[b]) - abs(Vi[b])
    sv = sig[nb-1+i]
    push!(out, (bus = b, name = name, quantity = :vm, multi_area = abs(V[b]), integrated = abs(Vi[b]), difference = dvm, sigma = sv === nothing ? NaN : sv, difference_sigma = sv === nothing ? NaN : dvm / sv))
    maxVm = max(maxVm, abs(dvm))
    if i != s
      dva = rad2deg(angle(V[b] / Vi[b]))
      sa = sig[vacol(i)]
      sad = sa === nothing ? NaN : rad2deg(sa)
      push!(out, (bus = b, name = name, quantity = :va, multi_area = rad2deg(angle(V[b])), integrated = rad2deg(angle(Vi[b])), difference = dva, sigma = sad, difference_sigma = sa === nothing ? NaN : dva / sad))
      maxVa = max(maxVa, abs(dva))
    end
  end
  rel = _se_tap_release_compare(tap_rows, res.tapEstimates)
  return (rows = out, max_dvm_pu = maxVm, max_dva_deg = maxVa, integrated_converged = res.converged, tap_release = rel, tap_release_differences = count(t -> !t.same, rel))
end

"""
    runse_multi_area!(net, measurements = net.measurements; areas = nothing,
                      partition = nothing, max_passes = nothing,
                      compare_integrated = nothing, update_net = true) -> MultiAreaSEResult

Multi-area state estimation, two-level method. The net is split into areas
([`area_partition`](@ref): from the map `areas` bus name => area, from an
area list file ([`read_area_list`](@ref)) when `areas` is a path, else
from `state_estimation.multi_area.partition`, `bus_area` (the stored bus
areas) by default; or pass a ready `partition`). `max_passes` and
`compare_integrated` default to the keys of `state_estimation.multi_area`. Every row belongs to one
area: a bus row to the area of its bus, a branch row to the area of its
measured end. Per pass each area estimates its interior with its own rows
and the border buses held at the coordinator's values; the coordinator then
takes one Gauss-Newton step on the border states from the sum of the area
models (Schur complement of each area's gain matrix onto its border columns
and its reduced gradient). The fixed point is the integrated estimate on the
same rows; the passes stop when no border state moves by more than
`state_estimation.tol`.

An area whose own rows do not determine its interior for given border
values is refused before the first pass, with its rank deficits. Not
supported (refused with the reason): more than one slack bus,
`state_estimation.zero_injection = constraint`, PMU angle rows in more than
one area.

With `compare_integrated` the result carries the integrated estimate on the
same rows and the difference per state, also in units of its standard
deviation. With `update_net` the estimated voltages are written to the
buses of `net`.
"""
function runse_multi_area!(net::Net, measurements::Vector{Measurement} = Measurement[m for m in net.measurements]; areas::Union{Nothing,AbstractDict,AbstractString} = nothing, partition::Union{Nothing,AreaPartition} = nothing, max_passes::Union{Nothing,Integer} = nothing, compare_integrated::Union{Nothing,Bool} = nothing, update_net::Bool = true)::MultiAreaSEResult
  cfg = state_estimation_config()
  part = partition === nothing ? _se_configured_partition(net, cfg, areas) : partition
  max_passes = something(max_passes, cfg.multi_area.max_passes)
  compare_integrated = something(compare_integrated, cfg.multi_area.compare_integrated)
  # the range check of the voltage rows, as in runse!
  checked, rc = _se_range_check(net, measurements)
  _se_log_range_check(net, rc, "runse_multi_area!")
  ma = _runse_two_level(net, checked, cfg, part; max_passes = Int(max_passes), check = true)
  return _multi_area_result(net, checked, cfg, part, ma; compare_integrated = compare_integrated, update_net = update_net)
end

function _multi_area_result(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig, part::AreaPartition, ma; compare_integrated::Bool, update_net::Bool)
  V = ma.V
  jsum, dof = _se_multi_area_j_dof(part, ma)
  arows = NamedTuple[]
  borderset = Set(part.border_buses)
  for a in part.areas
    run = ma.runs[a]
    ja = sum(run.model.w .* run.model.r .^ 2)
    push!(arows, (area = a, name = part.area_names[searchsortedfirst(part.areas, a)], buses = run.sub.own, border_buses = [b for b in run.sub.own if b in borderset], rows = length(ma.rows[a]), j = ja, iterations = run.result.iterations, converged = run.result.converged, result = run.result))
  end
  brows = NamedTuple[(bus = b, name = getCompName(net.nodeVec[b].comp), area = part.area_of_bus[b], vm_pu = abs(V[b]), va_deg = rad2deg(angle(V[b]))) for b in part.border_buses]
  cmp = compare_integrated ? _se_multi_area_compare(net, measurements, cfg, part, V; tap_rows = ma.tap_rows) : nothing
  if cmp !== nothing && cmp.tap_release_differences > 0
    diff = [string(t.name, " (area ", t.area, "): ", t.multi_area[1], "/", t.multi_area[2], ", integrated ", t.integrated[1], "/", t.integrated[2]) for t in cmp.tap_release if !t.same]
    @warn "runse_multi_area!: the tap release of $(length(diff)) regulator(s) differs from the integrated estimate: " * join(diff, "; ")
  end
  if update_net
    for (b, nd) in enumerate(net.nodeVec)
      nd._vm_pu = abs(V[b])
      nd._va_deg = rad2deg(angle(V[b]))
    end
  end
  return MultiAreaSEResult(V, ma.converged, ma.passes, ma.history, jsum, dof, arows, brows, part, isempty(ma.tap_rows) ? nothing : ma.tap_rows, cmp)
end

# ---------------------------------------------------------------------------
# Bad data (decision E6)
# ---------------------------------------------------------------------------

## Normalized residuals of the active rows of a converged two-level run with
## the exact residual covariance of the whole estimate, Omega = R - H C H'.
## The covariance of the border variables is the inverse of the coordinator
## matrix, C_BB = (sum_a S_a)^-1; for a row i of area a with interior part
## h_I and border part h_B (the area's border variables, the slack angle
## held fixed):
##   h' C h = h_I' G_II^-1 h_I + t' C_BB t,   t = h_B - G_BI G_II^-1 h_I,
## which follows from C_II = G_II^-1 + X C_BB X', C_IB = -X C_BB,
## X = G_II^-1 G_IB. Keyed by the row's index in `measurements`; LINKAGG
## aggregates (no single source row) are left out, as in `runse_diagnostics`.
function _se_multi_area_residuals(ma)
  nv = length(ma.vars)
  areas = sort!(collect(keys(ma.runs)))
  A, _ = _se_coordinator_system([ma.runs[a].model for a in areas], ma.vix, nv)
  Cbb = inv(Symmetric(A))
  out = Dict{Int,NamedTuple}()
  for a in areas
    run = ma.runs[a]
    md = run.model
    idx = [get(ma.vix, v, 0) for v in md.vars]
    nb = length(idx)
    Ca = zeros(nb, nb)
    for p = 1:nb, q = 1:nb
      (idx[p] > 0 && idx[q] > 0) && (Ca[p, q] = Cbb[idx[p], idx[q]])
    end
    Z = md.F === nothing ? nothing : md.F \ Matrix(md.HI')
    for i in eachindex(md.w)
      wi = md.w[i]
      wi > 0.0 || continue
      pos = run.lin.orig[i]
      pos == 0 && continue
      hB = md.HB[i, :]
      q = if Z === nothing
        dot(hB, Ca * hB)
      else
        zi = Z[:, i]
        t = hB .- md.GIB' * zi
        dot(Vector(md.HI[i, :]), zi) + dot(t, Ca * t)
      end
      om = max(1.0 / wi - q, 0.0)
      rn = om <= eps(Float64) ? 0.0 : md.r[i] / sqrt(om)
      out[run.kept[pos]] = (normalized_residual = rn, wii = om * wi, residual = md.r[i])
    end
  end
  return out
end

"""
    runse_multi_area_diagnostics(net, measurements = net.measurements; areas = nothing,
                                 partition = nothing, max_passes = nothing) -> NamedTuple

Bad-data elimination for the multi-area estimation, by the rule of
[`runse_diagnostics`](@ref): while J fails the band test above, the row with
the largest normalized residual (at least `state_estimation.k_eliminate`,
localizable, not protected from elimination, not a range-check substitute) is
switched off and the two-level estimation runs again, one row per round, up
to `state_estimation.max_eliminations`. The normalized residuals use the
exact residual covariance of the whole estimate, the border information of
every area included, so a row at an area border is judged as in the
integrated estimate (decision E6: the largest over all areas wins, its owner
removes it).

Returns `(result, trace, stop_reason, measurements)`: the
[`MultiAreaSEResult`](@ref) of the last run, one trace row per elimination
(`elimination`, `measurement_index`, `id`, `typ`, `area`,
`normalized_residual_before`, `residual_before`, `wii`, `objective_before`,
`objective_after`), the stop reason (`:consistent`, `:band_low`,
`:no_suspicious_left`, `:no_localizable_suspect`, `:max_eliminations`,
`:not_converged`) and the measurement vector with the eliminated rows
inactive.
"""
function runse_multi_area_diagnostics(net::Net, measurements::Vector{Measurement} = Measurement[m for m in net.measurements]; areas::Union{Nothing,AbstractDict,AbstractString} = nothing, partition::Union{Nothing,AreaPartition} = nothing, max_passes::Union{Nothing,Integer} = nothing)
  cfg = state_estimation_config()
  part = partition === nothing ? _se_configured_partition(net, cfg, areas) : partition
  max_passes = something(max_passes, cfg.multi_area.max_passes)
  checked, rc = _se_range_check(net, measurements)
  _se_log_range_check(net, rc, "runse_multi_area_diagnostics")
  substituted = Set{Int}(r.measurement_index for r in rc if r.reason === :range_check_substituted)
  meas2 = copy(checked)
  trace = NamedTuple[]
  stop = :consistent
  ma = nothing
  firstRound = true
  while true
    sparlectra_check_abort()
    ma = _runse_two_level(net, meas2, cfg, part; max_passes = Int(max_passes), check = firstRound)
    firstRound = false
    j, dof = _se_multi_area_j_dof(part, ma)
    isempty(trace) || (trace[end] = merge(trace[end], (objective_after = j,)))
    if !ma.converged
      stop = :not_converged
      break
    end
    verdict = _band_test_verdict(j, dof)
    if verdict.passed
      stop = :consistent
      break
    end
    if verdict.reason === :low
      stop = :band_low
      break
    end
    res = _se_multi_area_residuals(ma)
    ranked = sort!([(k, v) for (k, v) in res]; by = kv -> -abs(kv[2].normalized_residual))
    suspects = [(k, v) for (k, v) in ranked if abs(v.normalized_residual) >= cfg.k_eliminate && !_elimination_protected(meas2[k]) && !(k in substituted)]
    candidates = [(k, v) for (k, v) in suspects if v.wii > _SE_WII_THRESHOLD]
    if isempty(candidates)
      stop = isempty(suspects) ? :no_suspicious_left : :no_localizable_suspect
      break
    end
    if length(trace) >= cfg.max_eliminations
      stop = :max_eliminations
      break
    end
    k, v = candidates[1]
    m = meas2[k]
    push!(trace, (elimination = length(trace) + 1, measurement_index = k, id = m.id, typ = m.typ, area = _se_area_owner(net, part, m), normalized_residual_before = v.normalized_residual, residual_before = v.residual, wii = v.wii, objective_before = j, objective_after = NaN))
    meas2[k] = _set_measurement_active(m, false)
  end
  result = _multi_area_result(net, meas2, cfg, part, ma; compare_integrated = false, update_net = false)
  return (result = result, trace = trace, stop_reason = stop, measurements = meas2)
end

## J over all rows and the degrees of freedom of a two-level run (rows of
## the areas disjoint; states: the interiors plus the border variables).
function _se_multi_area_j_dof(part::AreaPartition, ma)
  j = 0.0
  m = 0
  nint = 0
  for a in part.areas
    md = ma.runs[a].model
    j += sum(md.w .* md.r .^ 2)
    m += count(>(0.0), md.w)
    nint += md.ninterior
  end
  return j, m - nint - length(ma.vars)
end
