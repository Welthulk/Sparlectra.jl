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

# file: src/api/run_state_estimation_service.jl
# purpose: Web UI/service state-estimation run (SE phase 5): build the net
#          through the shared import paths, read a measurement CSV v1 set,
#          check observability, run runse_diagnostics plus the final
#          runse!(updateNet = true), and persist the SE artifacts including
#          the se_state.csv chain anchor. Also the SE-started power-flow
#          chain run (_run_pf_from_se_service). Mirrors
#          run_contingency_service.jl (mode flag on POST /powerflow/run,
#          artifacts collected from the run dir).

"""
    _se_import_case(case_path, config) -> ImportedCase

Import a case for the state-estimation workflow through the one import
entry point. Throws on unsupported formats and on import errors. Shared by
the SE service run and the Web UI measurement generator so both accept
exactly the same cases, and both run with `ImportedCase.config` (the
CGMES start-value decision and the auto-profile rewrites reach the SE
solves exactly like the power-flow service, step 3a of the adapter task).
"""
function _se_import_case(case_path::AbstractString, config; requested_format::Symbol = :auto)::ImportedCase
  format = _detect_case_format(String(case_path); requested = requested_format)
  # format policy stays here (the wording is this service's contract)
  format in (:scf, :matpower, :cgmes, :dtf_for001) || throw(ArgumentError("State estimation needs a MATPOWER, CGMES, Sparlectra Case Format or DTF case; got format $(format)."))
  return import_case(String(case_path), config; requested_format = requested_format, run_kind = :state_estimation)
end

_se_import_case_net(case_path::AbstractString, config; requested_format::Symbol = :auto)::Net = _se_import_case(case_path, config; requested_format).net

## Format for decisions that only ask "which kind of case is this" (bus
## reference in the measurement CSV, CGMES specifics). A bare `.DAT` is
## ambiguous without an explicit request (FOR001 network vs FOR002
## reference), and that ambiguity must not abort a run that does not even
## depend on the answer, so it resolves to :unknown here.
function _se_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  return try
    _detect_case_format(String(case_path); requested = requested)
  catch
    :unknown
  end
end

## A phase tap changer counts as DECLARED nameplate data (not a constructor
## default) when the import stored its regulating-vector direction psi in
## tap_est_alpha_deg (mpc.sparlectra.tap_changers, optional psi column,
## default 90). Only declared phase changers join the mass release and the
## generator's deviation targeting.
_declared_phase_tap(br)::Bool = br.has_phase_tap && (br.phase_step_deg > 0.0 || br.phase_du_step > 0.0) && br.tap_est_alpha_deg != 0.0

## from-run truth resolution (measurement generator v2): resolve `run_id` in
## the persistent run index under `run_root`, enforce the case binding, and
## load the solved voltages into `net` WITHOUT re-solving. SE runs read
## se_state.csv (bus by name); PF runs read the detail CSV
## bus_voltages_complex.csv, which exists only when the run wrote detail
## CSVs. A missing artifact rejects up front with the hint to repeat the run
## with the detail CSV enabled; there is deliberately no silent re-solve.
function _se_truth_from_run!(net::Net, run_root::AbstractString, run_id::AbstractString, case_path::AbstractString)
  index = load_powerflow_run_index(run_root)
  entry = nothing
  for e in get(index, "runs", Any[])
    e isa AbstractDict || continue
    if String(get(e, "run_id", "")) == run_id
      entry = e
      break
    end
  end
  entry === nothing && throw(ArgumentError("run $(run_id) not found in the run history"))
  get(entry, "success", false) == true || throw(ArgumentError("run $(run_id) was not successful; the generator only accepts converged runs"))
  bound = basename(String(get(entry, "casefile", "")))
  bound == basename(case_path) || throw(ArgumentError("run $(run_id) belongs to case $(isempty(bound) ? "?" : bound), not $(basename(case_path))"))
  run_mode = String(get(entry, "run_mode", ""))
  run_mode in ("", "se") || throw(ArgumentError("run $(run_id) is a $(run_mode) run; only power-flow and state-estimation runs provide a truth state"))
  outdir = String(get(entry, "output_dir", joinpath(run_root, run_id)))
  ts = String(get(entry, "timestamp", ""))
  if run_mode == "se"
    f = joinpath(outdir, "se_state.csv")
    isfile(f) || throw(ArgumentError("run $(run_id) has no se_state.csv artifact; repeat the state estimation"))
    readSEStateCSV!(net; file = f)
    return (kind = "se", timestamp = ts, source_file = "se_state.csv")
  end
  f = joinpath(outdir, "bus_voltages_complex.csv")
  isfile(f) || throw(ArgumentError("run $(run_id) has no bus_voltages_complex.csv artifact; repeat the run with the detailed result CSV enabled"))
  _se_apply_pf_voltage_csv!(net, f)
  return (kind = "pf", timestamp = ts, source_file = "bus_voltages_complex.csv")
end

## tolerant reader for the PF detail CSV bus_voltages_complex.csv: the file
## is written in a user-selected format (technical/excel_us: ',' delimiter,
## '.' decimals; excel_de: ';' delimiter, ',' decimals). Buses match by
## bus_name first, then by the original bus number; any shape or matching
## failure aborts with a line-precise error instead of guessing.
function _se_apply_pf_voltage_csv!(net::Net, file::AbstractString)
  lines = readlines(file)
  header_ln = findfirst(l -> !startswith(strip(l), "#") && !isempty(strip(l)), lines)
  header_ln === nothing && throw(ArgumentError("$(file): empty file"))
  delim = occursin(";", lines[header_ln]) ? ";" : ","
  cols = [strip(String(c)) for c in split(lines[header_ln], delim)]
  bi = findfirst(==("bus"), cols)
  ni = findfirst(==("bus_name"), cols)
  vi = findfirst(==("vm_pu"), cols)
  ai = findfirst(==("va_deg"), cols)
  (bi === nothing || vi === nothing || ai === nothing) && throw(ArgumentError("$(file): missing bus/vm_pu/va_deg columns"))
  num(s) = delim == ";" ? tryparse(Float64, replace(strip(s), "," => ".")) : tryparse(Float64, strip(s))
  orig_to_idx = Dict{Int,Int}(o => i for (i, o) in net.busOrigIdxDict)
  n = length(net.nodeVec)
  seen = falses(n)
  for ln = (header_ln+1):length(lines)
    line = strip(lines[ln])
    isempty(line) && continue
    fields = split(line, delim)
    length(fields) == length(cols) || throw(ArgumentError("$(file):$(ln): expected $(length(cols)) fields, got $(length(fields)) (a name containing the delimiter?)"))
    idx = nothing
    if ni !== nothing
      name = String(strip(fields[ni]))
      haskey(net.busDict, name) && (idx = net.busDict[name])
    end
    if idx === nothing
      bno = tryparse(Int, strip(fields[bi]))
      bno !== nothing && (idx = get(orig_to_idx, bno, (1 <= bno <= n && isempty(orig_to_idx)) ? bno : nothing))
    end
    idx === nothing && throw(ArgumentError("$(file):$(ln): bus '$(strip(fields[bi]))' not found in the case"))
    vm = num(String(fields[vi]))
    va = num(String(fields[ai]))
    (vm === nothing || va === nothing) && throw(ArgumentError("$(file):$(ln): non-numeric vm_pu/va_deg"))
    setVmVa!(node = net.nodeVec[idx], vm_pu = vm, va_deg = va)
    seen[idx] = true
  end
  for i = 1:n
    seen[i] || getNodeType(net.nodeVec[i]) == Isolated || throw(ArgumentError("$(file): no voltage row for bus $(i); the file does not match this case"))
  end
  return nothing
end

"""
    _seeded_permutation(rng, v) -> Vector

A random permutation of `v` that depends only on the RNG's Float64 stream.
`Random.shuffle` changed its algorithm in Julia 1.13: the same
`MersenneTwister` seed then drew other rows than on 1.12, so a generated
set with "bad data on 1 row" hit a different measurement per Julia version
(found 2026-09-11 when the extended profile went red on 1.13 only). The
Float64 stream of the Mersenne Twister is the same on both, so ordering
by one uniform draw per element keeps "the seed decides which rows" true
across Julia versions.
"""
_seeded_permutation(rng::Random.AbstractRNG, v::AbstractVector) = v[sortperm(rand(rng, length(v)))]

"""
    _se_generate_measurement_set(case_path, out_path; kwargs...) -> NamedTuple

Service backend of the Web UI measurement generator (the Web UI layer never
calls a solver directly): import the case through the shared SE import
paths, optionally shift up to `tap_count` seed-randomly selected
in-service transformers by `tap_steps` MECHANICAL steps each on their own
grids (whole steps only, a tap changer has no half positions; the Web UI
enforces this before the call; the seed decides WHICH transformers,
non-machine declared changers are drawn first), obtain the truth state
(`truth_source = :fresh_solve` solves the
configured power flow at reference tightness; `:from_run` adopts the
solved voltages of run `run_id` under `run_root` without re-solving and
is mutually exclusive with a tap deviation), optionally keep only one
balance-aware flow end per branch (`flow_ends`), rewrite passive-node
balances at `passive_sigma` or as protected zero-injection constraints
(`passive_as_zi`), generate a measurement CSV v1 with percent-of-reading
sigmas, optional seeded noise, optional currents/current angles, optional
gross errors (`gross_k` sigma on `gross_count` seed-randomly selected
telemetry rows; protected zero-injection rows are never corrupted, and
the seed decides WHICH rows), readable name-based ids, the structured tap
comment table, and the
in-file case binding. Returns `(rows, noisy, tap_note, gross_note,
island_note, truth_note, flow_note, passive_note)` for the caller's
confirmation message. Throws
`ArgumentError` with a user-readable text when the case has no in-service
transformer for a requested deviation or the power flow does not converge.
"""
function _se_generate_measurement_set(
  case_path::AbstractString,
  out_path::AbstractString;
  noise::Bool,
  gross_k::Float64,
  gross_count::Int = 1,
  tap_steps::Float64,
  tap_count::Int = 1,
  include_i::Bool,
  sigma_u_pct::Float64,
  sigma_i_pct::Float64,
  sigma_p_pct::Float64,
  sigma_q_pct::Float64,
  sigma_ia_deg::Float64,
  truth_source::Symbol = :fresh_solve,
  run_id::AbstractString = "",
  run_root::Union{Nothing,AbstractString} = nothing,
  flow_ends::Symbol = :both,
  passive_sigma::Float64 = 0.05,
  passive_as_zi::Bool = false,
  seed::Int = 42,
)
  truth_source in (:fresh_solve, :from_run) || throw(ArgumentError("truth source must be fresh_solve or from_run"))
  flow_ends in (:both, :one_balance_aware) || throw(ArgumentError("flow measurements per branch must be both or one_balance_aware"))
  passive_sigma > 0.0 || throw(ArgumentError("passive-node balance sigma must be positive"))
  gross_count >= 1 || throw(ArgumentError("bad data count must be at least 1"))
  tap_count >= 1 || throw(ArgumentError("tap deviation transformer count must be at least 1"))
  if truth_source === :from_run
    # a run state is a finished snapshot: the generator never re-solves it,
    # so a tap deviation (which needs a fresh solve of the shifted model)
    # is rejected here even if the GUI lock was bypassed
    tap_steps == 0.0 || throw(ArgumentError("tap deviation requires truth state 'fresh solve'; a run state is a finished snapshot and is not re-solved"))
    isempty(strip(run_id)) && throw(ArgumentError("truth state 'from run' needs a run id"))
    run_root isa AbstractString || throw(ArgumentError("truth state 'from run': no run directory root available"))
  end
  # the generator honors the case's own settings the same way a run does
  config = resolve_config(DEFAULT_SPARLECTRA_CONFIG_PATH, case_path).config
  imported = _se_import_case(case_path, config)
  net = imported.net
  config = imported.config
  # A generated set REPLACES; it is never an addition. A Sparlectra Case
  # Format case brings its own measurements along, and without this the
  # generator wrote them out again next to the fresh ones (59 carried + 75
  # generated = 134 rows in a file that should have had 75).
  empty!(net.measurements)
  tap_note = ""
  tap_branches = Int[]
  if tap_steps != 0.0
    # never PREFER a machine (generator step-up) transformer: the tap mass
    # release skips those by design, so a deviation there could never be
    # resolved by "estimate taps" (seen on a case where the FIRST trafo was
    # the GSU behind the slack: 320 MVar of unexplainable imbalance).
    # Selection: up to tap_count transformers, drawn SEED-randomly (the
    # user sets how many at most, the seed decides which) with the class
    # priority preserved: non-machine declared changers (ratio tap, or
    # nameplate phase tap of a PST) fill first, then declared changers on
    # machine trafos, then any transformer at all (the run warns there)
    pool1 = [k for k in eachindex(net.branchVec) if (br = net.branchVec[k]; br.ratio != 0.0 && (br.has_ratio_tap || _declared_phase_tap(br)) && _branch_terminal_state(br) == :closed && !_is_machine_transformer(net, k))]
    pool2 = [k for k in eachindex(net.branchVec) if (br = net.branchVec[k]; br.ratio != 0.0 && (br.has_ratio_tap || _declared_phase_tap(br)) && _branch_terminal_state(br) == :closed && !(k in pool1))]
    pool3 = [k for k in eachindex(net.branchVec) if (br = net.branchVec[k]; br.ratio != 0.0 && _branch_terminal_state(br) == :closed && !(k in pool1) && !(k in pool2))]
    isempty(pool1) && isempty(pool2) && isempty(pool3) && throw(ArgumentError("tap deviation requested but the case has no in-service transformer"))
    # a dedicated RNG stream keeps the transformer draw independent of the
    # noise draw: the same seed always hits the same transformers, even
    # when noise or sigma settings change between regenerations.
    # Draw from ONE class only, the best non-empty one: mixing estimable
    # and machine transformers would poison the set, because a machine
    # (generator step-up) deviation is skipped by the mass release BY
    # DESIGN and leaves an unexplainable model error (seen: J ~ 318
    # instead of ~ dof when max 2 spilled onto the GSU of a 2-trafo case).
    # tap_count is therefore capped at the chosen class; machine or
    # arbitrary transformers serve only cases with nothing estimable at
    # all, and the estimator run warns there.
    rng_tap = Random.MersenneTwister(seed * 4093 + 11)
    pool = !isempty(pool1) ? pool1 : (!isempty(pool2) ? pool2 : pool3)
    picks = _seeded_permutation(rng_tap, pool)
    tap_branches = sort!(picks[1:min(tap_count, length(picks))])
    tap_limit_note = length(tap_branches) < tap_count ? " (limited to $(length(tap_branches)) eligible transformer(s), max $(tap_count) requested)" : ""
    nby = _bus_name_by_idx(net)
    notes = String[]
    for ti in tap_branches
      br = net.branchVec[ti]
      if br.has_ratio_tap || !_declared_phase_tap(br)
        base_tap = br.tap_ratio == 0.0 ? br.ratio : br.tap_ratio
        # cascade convention (reciprocal from side): n steps on the fraction
        # grid mean tap_ratio = base / (1 + n * tap_step), so an integer n is
        # exactly a mechanical position of the estimator's fixation grid
        dev_step = br.tap_step > 0.0 ? br.tap_step : 0.00625
        br.tap_ratio = base_tap / (1.0 + tap_steps * dev_step)
        pct_note = " (~$(round((1.0 / (1.0 + tap_steps * dev_step) - 1.0) * 100.0; digits = 2))% ratio)"
        push!(notes, "tap deviation $(tap_steps) steps$(pct_note) on transformer branch $(ti) ($(get(nby, Int(br.fromBus), string(Int(br.fromBus))))-$(get(nby, Int(br.toBus), string(Int(br.toBus)))))")
      else
        # pure phase shifter (PST): the deviation lives on the changer's OWN
        # mechanical grid. Delta-u PST: n additional-voltage steps through
        # the cascade (the shift angle follows via atan, it is NOT the
        # grid); degree-grid PST: additive degrees.
        if br.phase_du_step > 0.0
          ψ = deg2rad(br.tap_est_alpha_deg)
          tbase = (br.ratio == 0.0 ? 1.0 : br.ratio) * cis(deg2rad(br.angle))
          tlive = (br.tap_ratio == 0.0 ? abs(tbase) : br.tap_ratio) * cis(deg2rad(br.phase_shift_deg))
          r2live = real((tbase / tlive - 1.0) * cis(-ψ))
          tnew = tbase / (1.0 + (r2live + tap_steps * br.phase_du_step) * cis(ψ))
          br.tap_ratio = abs(tnew)
          br.phase_shift_deg = rad2deg(angle(tnew))
          push!(notes, "phase deviation $(tap_steps) Delta-u step(s) ($(br.phase_du_step) pu each, psi $(br.tap_est_alpha_deg) deg) on PST branch $(ti) ($(get(nby, Int(br.fromBus), string(Int(br.fromBus))))-$(get(nby, Int(br.toBus), string(Int(br.toBus)))))")
        else
          br.phase_shift_deg += tap_steps * br.phase_step_deg
          push!(notes, "phase deviation $(tap_steps) step(s) ($(br.phase_step_deg) deg each) on PST branch $(ti) ($(get(nby, Int(br.fromBus), string(Int(br.fromBus))))-$(get(nby, Int(br.toBus), string(Int(br.toBus)))))")
        end
      end
    end
    tap_note = string(", ", join(notes, "; "), tap_limit_note)
  end
  truth_comment = ""
  truth_note = ""
  if truth_source === :fresh_solve
    # config-driven solve: picks up island-wise solving (multi-island CGMES
    # deliveries) and the other configured solver options. The tolerance is
    # tightened to 1e-8: the measurement values ARE the reference truth, and
    # a loosely converged state costs the SE-started chain PF its warm start.
    pf_cfg = _copy_powerflow_with(config.powerflow; tol = min(config.powerflow.tol, 1e-8), max_iter = max(config.powerflow.max_iter, 40))
    ite, erg = runpf!(net; config = pf_cfg)
    erg == 0 || throw(ArgumentError("power flow on $(basename(case_path)) did not converge; no measurements generated"))
    truth_comment = "truth: fresh solve, $(ite) iteration(s), tol $(pf_cfg.tol), island-wise"
    truth_note = ", pre-solve $(ite) iteration(s) at tol $(pf_cfg.tol)"
  else
    info = _se_truth_from_run!(net, String(run_root), String(strip(run_id)), case_path)
    truth_comment = "truth: run $(strip(run_id)) ($(info.kind), $(info.timestamp), $(info.source_file)); state adopted, not re-solved"
    truth_note = ", truth from $(info.kind) run $(strip(run_id))"
  end
  # measurements cover every island: the estimator solves island-wise
  # (per-island reference), so the whole delivery is estimable. Count on
  # the ELECTRICAL view (closed links connect): the raw net is not
  # link-contracted here, and the branch-only solver view would report a
  # busbar section joined by a closed coupler as its own ref-less island
  n_islands = length(electricalIslandComponents(net))
  island_note = n_islands > 1 ? ", $(n_islands) islands" : ""
  # seeded noise keeps the demo set reproducible across regenerations;
  # relativeSigma: the fractions below are percent-of-reading / 100.
  stddev = measurementStdDevs(vm = sigma_u_pct / 100.0, pinj = sigma_p_pct / 100.0, qinj = sigma_q_pct / 100.0, pflow = sigma_p_pct / 100.0, qflow = sigma_q_pct / 100.0, imag = sigma_i_pct / 100.0, ia = sigma_ia_deg > 0.0 ? sigma_ia_deg : 0.05)
  meas = generateMeasurementsFromPF(net; noise = noise, stddev = stddev, relativeSigma = true, includeImag = include_i, includeIa = sigma_ia_deg > 0.0, rng = Random.MersenneTwister(seed))
  # readable ids: the defaults carry internal indices ("Vm_bus_1"), which
  # diverge from the bus/branch names on CGMES nets and read like a bug
  # next to mRID location columns (WebUI sets only, the API keeps its
  # stable default scheme)
  idnby = _bus_name_by_idx(net)
  short(t) = replace(string(t), "Meas" => "")
  for (gi2, m) in enumerate(meas)
    newid = if m.busIdx !== nothing
      string(short(m.typ), "_", get(idnby, m.busIdx, string(m.busIdx)))
    elseif m.branchIdx !== nothing
      string(short(m.typ), "_", getCompName(net.branchVec[m.branchIdx].comp), "_", m.direction)
    else
      m.id
    end
    meas[gi2] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = newid, linkIdx = m.linkIdx)
  end
  # flow measurements per branch: :both keeps the full from/to pairs;
  # :one_balance_aware keeps exactly one flow-measurement group per branch.
  # The end choice is deterministic: prefer the end whose bus carries a
  # TELEMETRY injection measurement (ZI pseudo-rows do not count); when
  # both or neither do, the from end wins. Choices land in the set comments.
  flow_note = ""
  flow_comments = String[]
  if flow_ends === :one_balance_aware
    injBuses = Set{Int}(m.busIdx for m in meas if m.active && m.typ in (PinjMeas, QinjMeas) && m.busIdx !== nothing && !startswith(m.id, "ZI"))
    flowBranches = sort!(unique(Int[m.branchIdx for m in meas if m.branchIdx !== nothing && m.typ in (PflowMeas, QflowMeas, ImagMeas, IaMeas)]))
    chosenByBranch = Dict{Int,Symbol}()
    for k in flowBranches
      br = net.branchVec[k]
      fromInj = Int(br.fromBus) in injBuses
      toInj = Int(br.toBus) in injBuses
      chosen = fromInj == toInj ? :from : (fromInj ? :from : :to)
      chosenByBranch[k] = chosen
      reason = fromInj == toInj ? (fromInj ? "both ends have injection telemetry" : "no end has injection telemetry") : "injection telemetry at the $(chosen) bus"
      push!(flow_comments, string("flow_end,", k, ",", getCompName(br.comp), ",", chosen, ",", reason))
    end
    meas = [m for m in meas if !(m.branchIdx !== nothing && m.typ in (PflowMeas, QflowMeas, ImagMeas, IaMeas) && m.direction !== chosenByBranch[m.branchIdx])]
    flow_note = ", one flow end per branch (balance-aware)"
  end

  # passive nodes: buses without generation, load, and shunt either get
  # their Pinj/Qinj balance rows at the CONFIGURED sigma with an exact
  # zero value, or (checkbox) protected zero-injection constraints via
  # addZeroInjectionMeasurements! and NO normal injection rows on top.
  passive_note = ""
  passive_comment = ""
  passiveBuses = findPassiveBuses(net)
  if !isempty(passiveBuses)
    nby_p = _bus_name_by_idx(net)
    pnames = join((get(nby_p, b, string(b)) for b in passiveBuses), " ")
    passiveSet = Set(passiveBuses)
    if passive_as_zi
      meas = [m for m in meas if !(m.busIdx !== nothing && m.busIdx in passiveSet && m.typ in (PinjMeas, QinjMeas))]
      # passive_sigma is the user-facing knob for BOTH passive modes: as a
      # balance row it is the reading's sigma, as a zero-injection row it is
      # how hard the constraint binds. It is in MW, so 0.001 means 1 kW.
      # Without this the constraint mode ignored the setting and used the
      # built-in floor, and there was no way to loosen it from the Web UI.
      zi_sigma = max(passive_sigma, ZERO_INJECTION_SIGMA)
      addZeroInjectionMeasurements!(meas; net = net, busIdxs = passiveBuses, sigma = zi_sigma)
      passive_note = ", $(length(passiveBuses)) passive node(s) as protected zero-injection constraints at sigma $(zi_sigma) MW"
      passive_comment = "passive: zero-injection constraints (ZI, sigma $(zi_sigma) MW, elimination-protected) at $(pnames)"
    else
      for (mi, m) in enumerate(meas)
        (m.busIdx !== nothing && m.busIdx in passiveSet && m.typ in (PinjMeas, QinjMeas)) || continue
        z = 0.0 + (noise ? randn(Random.MersenneTwister(seed * 1009 + 7 * m.busIdx + (m.typ == PinjMeas ? 0 : 1))) * passive_sigma : 0.0)
        meas[mi] = Measurement(typ = m.typ, value = z, sigma = passive_sigma, active = m.active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx)
      end
      passive_note = ", $(length(passiveBuses)) passive node(s) at balance sigma $(passive_sigma)"
      passive_comment = "passive: balance rows Pinj/Qinj = 0 with sigma $(passive_sigma) at $(pnames)"
    end
  end

  gross_note = ""
  if gross_k > 0.0
    # bad data on gross_count SEED-randomly drawn telemetry rows (the user
    # sets how many, the seed decides which); protected zero-injection
    # constraints are never corrupted (they are elimination-protected by
    # design, a gross error there could not be worked off). A dedicated RNG
    # stream keeps the row draw independent of the noise draw.
    eligible = [i for (i, m) in enumerate(meas) if m.active && !startswith(m.id, "ZI")]
    if !isempty(eligible)
      rng_gross = Random.MersenneTwister(seed * 7919 + 13)
      picks = _seeded_permutation(rng_gross, eligible)[1:min(gross_count, length(eligible))]
      sort!(picks)
      gross_ids = String[]
      for gi in picks
        m0 = meas[gi]
        meas[gi] = Measurement(typ = m0.typ, value = m0.value + gross_k * m0.sigma, sigma = m0.sigma, active = m0.active, busIdx = m0.busIdx, branchIdx = m0.branchIdx, direction = m0.direction, id = m0.id, linkIdx = m0.linkIdx)
        push!(gross_ids, m0.id)
      end
      gross_note = ", bad data +$(gross_k) sigma on $(length(gross_ids)) row(s): $(join(gross_ids, " "))"
    end
  end
  append!(net.measurements, meas)
  # record the tap positions the set was generated from as a structured
  # table (the estimator later has to explain any deviation between these
  # and the model): electrical_step = continuous position of the generation
  # state, fixed_step = nearest mechanical step, transferred_step = the
  # mechanical position an RTU would report (equal to fixed_step here)
  nby2 = _bus_name_by_idx(net)
  bn2(b) = get(nby2, Int(b), string(Int(b)))
  trmrids = _transformer_mrids(net)
  # the case binding lives IN the file (first comment after the version
  # line): renaming or re-uploading the set keeps the association, and the
  # SE service refuses a set bound to a different case up front
  # summary lines stay at the top; the bulky per-branch and per-row
  # blocks (flow_end choices, truth values) go BEHIND the taps table so
  # the comment-limited readers (set info, tap-deviation parse) always
  # reach the table on large cases
  # The noise state belongs IN the set. A noise-free set is an IDEAL one:
  # its J is 0 by construction, which looks like a perfect estimate and is
  # nothing of the sort. Recording it here is what lets the case file, the
  # state-estimation page and the run say so instead of leaving the reader
  # to guess from a suspiciously small J.
  noise_comment = noise ? "noise: gaussian (sigma U $(sigma_u_pct)%, P $(sigma_p_pct)%, Q $(sigma_q_pct)% of reading)" : "noise: none (ideal values, J is 0 by construction)"
  gencomments = String["generator: v2", "seed: $(seed)", noise_comment, truth_comment, "flow_ends: $(flow_ends)"]
  isempty(passive_comment) || push!(gencomments, passive_comment)
  tailcomments = copy(flow_comments)
  # noise-free truth values per row (the value BEFORE noise and gross
  # error): the SE run reads these comments back and writes the
  # measured/truth/estimated delta file se_deltas.csv. Passive and ZI
  # balance rows have truth 0 by construction.
  truthById = Dict{String,Float64}()
  tmeas = generateMeasurementsFromPF(net; noise = false, stddev = stddev, relativeSigma = true, includeImag = include_i, includeIa = sigma_ia_deg > 0.0)
  for tm in tmeas
    newid = if tm.busIdx !== nothing
      string(short(tm.typ), "_", get(idnby, tm.busIdx, string(tm.busIdx)))
    elseif tm.branchIdx !== nothing
      string(short(tm.typ), "_", getCompName(net.branchVec[tm.branchIdx].comp), "_", tm.direction)
    else
      tm.id
    end
    truthById[newid] = tm.value
  end
  passiveTruthSet = passive_as_zi ? Set{Int}() : Set(passiveBuses)
  for m in meas
    v = startswith(m.id, "ZI") ? 0.0 : get(truthById, m.id, nothing)
    (m.busIdx !== nothing && m.busIdx in passiveTruthSet && m.typ in (PinjMeas, QinjMeas)) && (v = 0.0)
    v === nothing && continue
    push!(tailcomments, string("truth_value,", m.id, ",", repr(v)))
  end
  tapcomments = vcat(String["case: $(basename(case_path))"], gencomments, String["sparlectra-taps v1", "branch,mrid,name,from_bus,to_bus,neutral_ratio,step,electrical_step,fixed_step,transferred_step,generation_deviation_steps"])
  ntr = 0
  for (k, br) in enumerate(net.branchVec)
    br.ratio != 0.0 || continue
    ntr += 1
    # electrical step on the trafo's own mechanical grid, the same grid the
    # estimator's fixation uses, so an injected integer deviation reads
    # back exactly: ratio taps on the FRACTION grid (r = base/current - 1
    # over tap_step), declared PSTs on the PHASE grid (degrees off the
    # neutral shift over phase_step_deg; the step column then carries
    # phase_step_deg)
    if br.has_ratio_tap || !_declared_phase_tap(br)
      cur = br.tap_ratio == 0.0 ? br.ratio : br.tap_ratio
      el = br.tap_step > 0.0 ? round((br.ratio / cur - 1.0) / br.tap_step; digits = 3) : NaN
      stepcol = br.tap_step
    elseif br.phase_du_step > 0.0
      # Delta-u PST: position in additional-voltage steps (the shift angle
      # is a consequence, not the grid)
      ψt = deg2rad(br.tap_est_alpha_deg)
      tbaset = (br.ratio == 0.0 ? 1.0 : br.ratio) * cis(deg2rad(br.angle))
      tlivet = (br.tap_ratio == 0.0 ? abs(tbaset) : br.tap_ratio) * cis(deg2rad(br.phase_shift_deg))
      el = round(real((tbaset / tlivet - 1.0) * cis(-ψt)) / br.phase_du_step; digits = 3)
      stepcol = br.phase_du_step
    else
      el = round((br.phase_shift_deg - br.angle) / br.phase_step_deg; digits = 3)
      stepcol = br.phase_step_deg
    end
    fx = isnan(el) ? "" : string(Int(round(el, RoundNearestTiesAway)))
    push!(tapcomments, string(k, ",", get(trmrids, k, ""), ",", getCompName(br.comp), ",", bn2(br.fromBus), ",", bn2(br.toBus), ",", br.ratio, ",", stepcol, ",", isnan(el) ? "" : el, ",", fx, ",", fx, ",", k in tap_branches ? tap_steps : 0.0))
  end
  # keep the in-file case binding (and the generator provenance) even when
  # the case has no transformers at all
  ntr == 0 && (tapcomments = vcat(String["case: $(basename(case_path))"], gencomments, String["case has no transformer branches (no tap steps to estimate)"]))
  append!(tapcomments, tailcomments)
  writeMeasurementsCSV(net; file = out_path, headerComments = tapcomments, busReference = _se_case_format(case_path) === :cgmes ? :mrid : :name)
  return (rows = length(net.measurements), noisy = noise, tap_note = tap_note, gross_note = gross_note, island_note = island_note, truth_note = truth_note, flow_note = flow_note, passive_note = passive_note)
end

## True when at least one transformer of `net` has its tap released as an
## estimation state.
## The one statement the tap fallback makes, shared by all three surfaces
## that must carry it: the Web UI result summary, the tap table (which the
## fallback REPLACES, because estimated-looking positions that are model
## values are worse than no table), and se_diagnostics.md. One constant so
## the three cannot drift apart; a run whose tap positions are model values
## is otherwise indistinguishable from a successful tap estimation.
const _SE_TAP_FALLBACK_NOTE = "tap estimation did NOT converge: every tap was frozen at its MODEL position, so the tap positions in this result are model values, not estimates. The J of this run therefore measures the model tap positions, and a large J is expected when those positions are wrong."

_any_tap_released(net::Net)::Bool = any(b -> b.tap_est_mode !== :none, net.branchVec)

## Freeze every released tap back to its model position and report how many
## were frozen. Used by the service fallback: a set too thin to pin its taps
## should still deliver an estimate of the voltages.
function _freeze_all_tap_estimation!(net::Net)::Int
  frozen = 0
  for (i, b) in enumerate(net.branchVec)
    b.tap_est_mode === :none && continue
    setTapEstimation!(net; trafo = i, enabled = false)
    frozen += 1
  end
  return frozen
end

## Iterations a run actually needed, not just its LAST solve. With released
## taps the run solves twice: the estimation with the tap states, then the
## fixation run with the taps nailed to their mechanical step. The second is
## short (it starts from a converged state), so reporting only that one hid
## the expensive half: a CGMES run needed 36 to 40 iterations first and
## reported "3", which is how an iteration cap of 30 could look sufficient
## while it broke that run (maintainer, 2026-09-06).
function _se_reported_iterations(res)::Int
  base = res.iterations
  tf = res.tapFixation
  tf === nothing && return base
  return max(base, Int(get(tf, :iterations_estimation, 0)), Int(get(tf, :iterations_fixation, 0)))
end

## shared import for the SE services (same paths as the other services);
## returns the full ImportedCase so the run continues on its effective
## config, and logs the CGMES start decision like the power-flow service
function _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata; requested_format::Symbol = :auto)
  format = _se_case_format(case_path; requested = requested_format)
  # named separately: "unknown" alone would leave the caller guessing, and
  # the way out (naming the format) is not obvious from the file name
  format === :unknown && return nothing, format, _api_failure("se_unsupported_format", "The case format could not be determined. A bare .DAT is ambiguous (FOR001 network case vs FOR002 reference file); pass case_format = :dtf_for001 for a DTF network case.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  format in (:matpower, :cgmes, :scf, :dtf_for001) || return nothing, format, _api_failure("se_unsupported_format", "State estimation needs a MATPOWER, CGMES, Sparlectra Case Format or DTF case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  imported = try
    _se_import_case(case_path, config; requested_format = requested_format)
  catch err
    err isa PowerFlowAborted && rethrow()
    return nothing, format, _api_failure("import_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  decision = get(imported.provenance, "cgmes_start_decision", nothing)
  decision === nothing || open(logfile, "a") do io
    println(io, decision)
  end
  return imported, format, nothing
end

"""
    _run_state_estimation_service(case_path, config_file, output_dir, run_id, measurement_file; kwargs...) -> SparlectraApiResult

Service backend of the Web UI "Run state estimation" action (SE phase 5).
Builds the net through the shared import paths, reads `measurement_file`
(measurement CSV v1, atomic), evaluates global observability (structural
islands plus FD-aware rank, phase 4), runs the bad-data diagnostics
(`runse_diagnostics`, sequential elimination on the configured budget) and
the final `runse!(updateNet = true)`, and writes the SE artifacts:
`measurements.csv` (copy), `se_diagnostics.md`, `se_view.md`,
`shunt_estimates.csv` (when shunts were released), and `se_state.csv` (the
chain anchor a later SE-started power flow consumes via `readSEStateCSV!`).

Options: `max_iter`, `tol`, `flatstart`, `robust`, `max_eliminations`,
`update_shunts`, `report_correlation` mirror the estimator keywords.
`tap_estimation = true` releases the tap of every in-service transformer
that carries a ratio tap changer (`setTapEstimation!` mode `:ratio`) before
the solve; the estimator then fixes each tap to its nearest mechanical step
and reports J before versus after the fixation (`se_tap_estimates.csv`).

Failure reasons: `se_unsupported_format`, `import_error`,
`invalid_measurements` (missing/unreadable/rejected file),
`se_not_observable` (observability quality `:not_observable`),
`se_not_converged`, plus the shared config failure.
"""
## headline mapping for the state-estimation timing file
const _SE_PERF_HEADLINE = (:case_loading_network_solver => "importing_case", :solver => "state_estimation", :postprocessing => "postprocessing_result", :artifact_writing => "writing_artifacts")

function _run_state_estimation_service(
  case_path::AbstractString,
  config_file::AbstractString,
  output_dir::AbstractString,
  run_id::String,
  measurement_file::AbstractString;
  # Solver settings default to `nothing`, NOT to a literal: a number here
  # would be a second source of truth next to the configuration, and it was
  # exactly that (service 6.0 against configuration 4.0 for k_suppress,
  # service 30 against configuration 20 for max_iter) which made the same
  # run behave differently depending on the entry point. `nothing` means
  # "not stated by the caller", and the run resolves it below against the
  # effective configuration, so the documented precedence holds: request,
  # then case configuration, then general configuration.
  max_iter::Union{Nothing,Int} = nothing,
  tol::Union{Nothing,Float64} = nothing,
  flatstart::Union{Nothing,Bool} = nothing,
  robust::Bool = false,
  max_eliminations::Union{Nothing,Int} = nothing,
  update_shunts::Bool = false,
  report_correlation::Bool = false,
  tap_estimation::Bool = false,
  k_eliminate::Union{Nothing,Float64} = nothing,
  robust_mode::Union{Nothing,Symbol} = nothing,
  robust_k1::Union{Nothing,Float64} = nothing,
  robust_k2::Union{Nothing,Float64} = nothing,
  k_suppress::Union{Nothing,Float64} = nothing,
  suppression_sigma::Union{Nothing,Float64} = nothing,
  case_format::Symbol = :auto,
  phase_callback = phase -> nothing,
)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "se")

  # Phase timing, same machinery and same file name the power-flow path uses.
  # Without it a slow estimation could only be guessed at: the topology
  # precheck once cost 45 s of a 55 s run on a 25000-bus network and nothing
  # in the run directory showed it. Phases are coarse on purpose (one line
  # per step, not per iteration), matching the Web UI instrumentation rule.
  phase_recorder = PowerFlowPhaseTimingRecorder()
  # The recorder alone only fills the result metadata. The Web UI status page
  # reads the JOB's phase, so every phase has to be announced as well - the
  # power-flow path does that through the same callback. Without it the page
  # froze at the last phase the service layer set (preparing_configuration).
  se_phase = function (name::AbstractString)
    _start_service_phase!(phase_recorder, name)
    try
      phase_callback(String(name))
    catch err
      # a reporting problem must never take the run down
      err isa PowerFlowAborted && rethrow()
      @debug "state estimation: phase callback failed" phase = name exception = err
    end
    return nothing
  end
  se_total_start = time_ns()

  # the same precedence the power-flow path uses (resolve_config, D5):
  # case configuration file, the case file's deprecated block, general
  # file, defaults
  se_phase("preparing_configuration")
  config = try
    resolve_config(config_file, case_path).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  # what the caller did not state comes from the EFFECTIVE configuration
  # (case file first, then the general one), never from a literal in this
  # signature
  se_cfg = config.state_estimation
  max_iter = something(max_iter, se_cfg.max_iter)
  tol = something(tol, se_cfg.tol)
  flatstart = something(flatstart, se_cfg.flatstart)
  k_eliminate = something(k_eliminate, se_cfg.k_eliminate)
  robust_mode = something(robust_mode, se_cfg.robust_mode)
  robust_k1 = something(robust_k1, se_cfg.robust_k1)
  robust_k2 = something(robust_k2, se_cfg.robust_k2)
  k_suppress = something(k_suppress, se_cfg.k_suppress)
  suppression_sigma = something(suppression_sigma, se_cfg.suppression_sigma)
  max_eliminations = something(max_eliminations, se_cfg.max_eliminations)

  # bad-data threshold surface (GUI-exposed): validated HERE, after the
  # values are resolved, with a user-readable failure instead of a deep
  # solver error. Before, this ran on the keywords alone, so a bad value in
  # a configuration file reached the solver unchecked.
  robust_mode in (:off, :staged, :replacement) || return _api_failure("invalid_request", "se_robust_mode must be off, staged, or replacement (got $(robust_mode)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  k_eliminate > 0.0 || return _api_failure("invalid_request", "se_k_eliminate must be positive.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  (robust_k1 > 0.0 && robust_k2 >= robust_k1) || return _api_failure("invalid_request", "se_robust_k1 must be positive and se_robust_k2 >= se_robust_k1 (got k1=$(robust_k1), k2=$(robust_k2)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  max_eliminations >= 0 || return _api_failure("invalid_request", "se_max_eliminations must be >= 0 (got $(max_eliminations)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  (k_suppress > 0.0 && suppression_sigma > 0.0) || return _api_failure("invalid_request", "se_k_suppress and se_suppression_sigma must be positive.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)


  se_phase("importing_case")
  imported, format, failure = _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata; requested_format = case_format)
  failure !== nothing && return failure
  net = imported.net
  config = imported.config

  # A Sparlectra Case Format case carries its measurements WITH the model
  # (#342), so a run on one needs no separate CSV: the set that came with the
  # case is used, and the run still writes measurements.csv as its artifact so
  # every SE run stays reproducible from its own output.
  from_case_file = format === :scf && !isfile(measurement_file) && !isempty(net.measurements)
  summary = if from_case_file
    open(logfile, "a") do io
      println(io, "measurements: ", length(net.measurements), " row(s) taken from the case file itself (no separate measurement set given)")
    end
    writeMeasurementsCSV(net; file = joinpath(output_dir, "measurements.csv"), busReference = :name)
    (total = length(net.measurements),)
  else
    if !isfile(measurement_file)
      hint = format === :scf ? " The case file carries no measurements either: generate a set, or export the case with its measurements." : ""
      return _api_failure("invalid_measurements", "Measurement file not found: $(measurement_file)." * hint, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    # case-binding gate: a set recorded for a different case would fail rows
    # deep in the reader ("bus X not found"); refuse it up front with the
    # actual cause. Sets without a binding (older files) pass with a log note.
    bound_case = _webui_measurement_set_case(String(measurement_file))
    # A case exported to SCF keeps the model of its source case, and it records
    # WHICH source it came from. A set bound to that source therefore fits this
    # run; without this the export silently invalidates every set the user
    # already generated.
    source_of_case = format === :scf ? _scf_source_reference(case_path) : ""
    if !isempty(bound_case) && bound_case != basename(String(case_path)) && bound_case != source_of_case
      return _api_failure("invalid_measurements", "Measurement set $(basename(String(measurement_file))) is bound to case $(bound_case) (recorded in the file), but this run uses $(basename(String(case_path))). Pick the set generated for this case, or regenerate one.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    isempty(bound_case) && open(logfile, "a") do io
      println(io, "note: the measurement set carries no case binding (# case: comment); older set, association checked by name only")
    end
    read_summary = try
      readMeasurementsCSV!(net; file = measurement_file, replace = true)
    catch err
      return _api_failure("invalid_measurements", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    cp(measurement_file, joinpath(output_dir, "measurements.csv"); force = true)
    read_summary
  end

  # Carry a duplicate finding into the run: it is the one property of a set that
  # explains a large J without any single measurement looking wrong.
  # Checked on the loaded set, not on the reader result: measurements carried
  # inside a case file take the other branch above and can be doubled just the
  # same.
  duplicate_rows = sum(values(_duplicate_measured_quantities(net.measurements)); init = 0)
  if duplicate_rows > 0
    base_metadata["se_duplicate_rows"] = duplicate_rows
    open(logfile, "a") do io
      println(io, "warning: $(duplicate_rows) of $(length(net.measurements)) measurements repeat a quantity that is already measured; expect an inflated J")
    end
  end

  # observability first (phase 4 two-stage verdict); a not-observable set is
  # a clean rejection, not a solver crash
  obs = try
    evaluate_global_observability(net)
  catch err
    return _api_failure("invalid_measurements", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  base_metadata["se_observability_quality"] = String(obs.quality)
  base_metadata["se_structural_islands"] = :structural_islands in obs.notes
  if hasproperty(obs, :islands)
    base_metadata["se_islands_total"] = length(obs.islands)
    base_metadata["se_islands_measured"] = obs.n_measured_islands
    base_metadata["se_islands_unmeasured"] = obs.n_unmeasured_islands
  end
  if obs.quality == :not_observable
    msg = if hasproperty(obs, :islands)
      bad = [string("island ", i.island, " (", i.n_bus, " buses): ", String(i.result.quality)) for i in obs.islands if i.measured && i.result.quality == :not_observable]
      "The measurement set does not observe every measured island ($(join(bad, "; "))). Each island is estimated with its own reference; add measurements to the failing island or deactivate its rows."
    else
      "The measurement set does not observe the system (quality :not_observable)."
    end
    return _api_failure("se_not_observable", msg, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # released transformer taps: every in-service transformer with a
  # ratio tap changer becomes an estimation state (mode :ratio; PST release
  # stays an API-level choice via setTapEstimation!). Released AFTER the
  # observability verdict: tap states are deliberately outside the
  # observability count, like the gated current rows.
  # the generator documents an injected tap deviation inside the set; a run
  # WITHOUT tap estimation then reports a large J that is a model
  # discrepancy, not bad data. Say so before anyone chases measurements.
  # both reads below mine the generator's COMMENTS in the measurement CSV
  # (documented tap deviation, truth values); a set that came with the case
  # file has no such file, so they are skipped rather than guessed
  # Measurements carried inside a case file must yield the same two facts a
  # CSV set yields, or an SCF run silently loses them: which taps the generator
  # deviated (the estimator has to absorb exactly those) and the per-row truth
  # values (the delta file).
  case_provenance = from_case_file ? _scf_measurement_provenance(case_path) : Dict{String,Any}()
  set_tap_devs = if from_case_file
    NamedTuple[(branch = Int(get(d, "branch", 0)), steps = Float64(get(d, "steps", 0.0))) for d in get(case_provenance, "tap_deviations", Any[]) if d isa AbstractDict]
  else
    _measurement_set_tap_deviations(String(measurement_file))
  end
  # A set that DOCUMENTS its tap deviations answers the question the
  # estimator would otherwise have to be told, so those transformers are
  # released automatically - exactly the ones named, nothing else. As a mere
  # hint this produced a J the user could not explain, and the hint appeared
  # only after a converged run. Bad-data elimination is the wrong tool here:
  # the measurements are correct and the model is not, so dropping rows
  # would hide a model error instead of fixing it.
  auto_tap_branches = Int[]
  if !tap_estimation && !isempty(set_tap_devs)
    for d in set_tap_devs
      b = d.branch
      (1 <= b <= length(net.branchVec)) || continue
      try
        setTapEstimation!(net; trafo = b, enabled = true)
        push!(auto_tap_branches, b)
      catch err
        open(logfile, "a") do io
          println(io, "note: tap estimation could not be released automatically on branch ", b, ": ", sprint(showerror, err))
        end
      end
    end
    open(logfile, "a") do io
      if isempty(auto_tap_branches)
        println(io, "note: the measurement set documents a tap deviation (", join(("$(d.steps) step(s) on branch $(d.branch)" for d in set_tap_devs), ", "), ") but no transformer could be released; the estimator cannot absorb it (expect band :high clustered at the transformer)")
      else
        println(io, "tap estimation released automatically on ", length(auto_tap_branches), " transformer(s) named by the measurement set: ", join(auto_tap_branches, ", "))
      end
    end
    base_metadata["se_set_tap_deviation"] = true
    base_metadata["se_auto_released_taps"] = auto_tap_branches
  end

  # topology precheck: the check RESULT is logged on every run,
  # the clean case included (a silent skip would read as a pass); findings
  # are advisory and never block the estimation
  se_phase("topology_precheck")
  topo_pre = validate_topology(net)
  # The findings are advisory and a large synthetic network produces many of
  # them (2170 on a 25000-bus case, 2155 of them "closed branch carries no
  # flow"). Printing every one buries the run log and nobody reads it, so the
  # log keeps a count per kind plus the first few of each; the complete list
  # stays in the metadata and in se_topology.csv.
  open(logfile, "a") do io
    println(io, topo_pre.summary)
    by_kind = Dict{Symbol,Vector{Any}}()
    for f in topo_pre.findings
      push!(get!(Vector{Any}, by_kind, f.kind), f)
    end
    for kind in sort!(collect(keys(by_kind)); by = string)
      fs = by_kind[kind]
      println(io, "  ", kind, ": ", length(fs), " finding(s)")
      for f in first(fs, 5)
        println(io, "    ", f.severity, " at ", f.location, ": ", f.evidence)
      end
      length(fs) > 5 && println(io, "    ... ", length(fs) - 5, " more (full list in se_topology.csv)")
    end
  end
  base_metadata["se_topology_findings"] = [Dict{String,Any}("stage" => String(f.stage), "kind" => String(f.kind), "location" => f.location, "evidence" => f.evidence, "severity" => String(f.severity)) for f in topo_pre.findings]
  base_metadata["se_topology_precheck_summary"] = topo_pre.summary

  tap_released = 0
  tap_skipped_machine = 0
  tap_machine_idxs = Int[]
  if tap_estimation
    for (k, br) in enumerate(net.branchVec)
      hasR = br.has_ratio_tap
      # DECLARED phase changers only (nameplate psi carried in
      # tap_est_alpha_deg by the mpc.sparlectra.tap_changers import): the
      # constructor gives every MATPOWER trafo default phase fields, and
      # mass-releasing those would add meaningless states everywhere
      hasP = _declared_phase_tap(br)
      (br.ratio != 0.0 && (hasR || hasP) && _branch_terminal_state(br) == :closed) || continue
      # machine (generator step-up) transformers are never mass-released:
      # the generator terminal voltage is set by the AVR, not observed, so
      # the tap state would absorb it. An explicit setTapEstimation! on a
      # specific transformer remains the deliberate way around this guard.
      if _is_machine_transformer(net, k)
        tap_skipped_machine += 1
        push!(tap_machine_idxs, k)
        continue
      end
      if hasR && hasP
        setTapEstimation!(net; trafo = k, mode = :both, alpha_deg = br.tap_est_alpha_deg)
      elseif hasP
        setTapEstimation!(net; trafo = k, mode = :pst, alpha_deg = br.tap_est_alpha_deg)
      else
        setTapEstimation!(net; trafo = k, mode = :ratio)
      end
      tap_released += 1
    end
    open(logfile, "a") do io
      tap_released == 0 && println(io, "tap estimation requested, but the case has no in-service transformer with a ratio or declared phase tap changer", tap_skipped_machine > 0 ? " (besides $(tap_skipped_machine) machine transformer(s), which are never mass-released)" : "")
      tap_skipped_machine > 0 && tap_released > 0 && println(io, "tap estimation: $(tap_skipped_machine) machine transformer(s) skipped (generator step-up; release explicitly via setTapEstimation! if intended)")
      # loud warning when the set's documented deviation sits on a SKIPPED
      # machine transformer: estimate taps can then never absorb it and J
      # stays large for a reason the user cannot see otherwise
      for d in set_tap_devs
        (1 <= d.branch <= length(net.branchVec) && _is_machine_transformer(net, d.branch)) || continue
        println(io, "WARNING: the set documents a tap deviation of $(d.steps) step(s) on branch $(d.branch), which is a machine transformer and SKIPPED by the mass release; the estimation cannot absorb it (release it explicitly via setTapEstimation!, or regenerate the set)")
      end
    end
  end

  # diagnostics (report) plus the final state-writing run (chain anchor)
  se_phase("estimation_diagnostics")
  diag = runse_diagnostics(net; max_eliminations = max_eliminations, maxIte = max_iter, tol = tol, flatstart = flatstart, robust = robust, reportResidualCorrelation = report_correlation, normalizedThreshold = k_eliminate, robustMode = robust_mode, robustK1 = robust_k1, robustK2 = robust_k2, kSuppress = k_suppress, suppressionSigma = suppression_sigma)
  # eliminated rows leave the FINAL run: the diagnostics identified them on
  # its own copy, so deactivate them here or the headline J would keep
  # carrying rows the workflow already removed (seen: an injected 10-sigma
  # error, eliminated with stop reason consistent, still pushed the
  # headline to J = 149 at dof 42 instead of J ~ dof)
  for e in diag.eliminations
    (1 <= e.measurement_index <= length(net.measurements)) || continue
    net.measurements[e.measurement_index] = _set_measurement_active(net.measurements[e.measurement_index], false)
  end
  se_phase("state_estimation")
  # kept for the tap fallback below: runse! writes its state into the net
  # whether or not it converged, so the retry needs the state as it was
  # BEFORE the first attempt, not the diverged iterate it gave up on
  v_before_taps = _any_tap_released(net) ? [(nd._vm_pu, nd._va_deg) for nd in net.nodeVec] : nothing
  res = runse!(net; maxIte = max_iter, tol = tol, flatstart = flatstart, updateNet = true, updateShunts = update_shunts, robust = robust, robustMode = robust_mode, robustK1 = robust_k1, robustK2 = robust_k2, kSuppress = k_suppress, suppressionSigma = suppression_sigma)

  # Released taps are extra states, and a measurement set that carries the
  # voltages fine can still be too thin to pin them: the estimate then does
  # not settle at all and the user gets nothing, although the SAME set
  # estimates cleanly without the taps (maintainer, 2026-09-06, case300 and
  # a CGMES delivery). So a non-convergence WITH released taps is not the
  # final answer: the taps are frozen back to their model position and the
  # estimation is repeated once. The log says it happened, because a silent
  # retry would hide that the reported taps are model values, not estimates.
  tap_fallback_used = false
  if !res.converged && _any_tap_released(net)
    frozen = _freeze_all_tap_estimation!(net)
    # The failed run wrote its state into the net (updateNet = true is not
    # gated on convergence), so without this the retry starts from the
    # DIVERGED iterate. Measured on the svedala set 2026-09-06: the failed
    # tap run left bus angles spread over -178 to +167 degrees where the
    # power flow had -26 to +40. With state_estimation.flatstart = true the
    # retry ignores the start state and the damage stays invisible, which is
    # why this never showed up; with a warm start it decides whether the
    # fallback converges at all.
    if v_before_taps !== nothing
      for (i, nd) in enumerate(net.nodeVec)
        nd._vm_pu, nd._va_deg = v_before_taps[i]
      end
    end
    open(logfile, "a") do io
      println(io, "state estimation did not converge with ", frozen, " released transformer tap(s); repeating WITHOUT tap estimation (the taps keep their model position). A set that cannot pin its taps needs more measurements around those transformers, not more iterations.")
    end
    @info "state estimation: tap estimation switched off after a non-converged run" released = frozen
    res = runse!(net; maxIte = max_iter, tol = tol, flatstart = flatstart, updateNet = true, updateShunts = update_shunts, robust = robust, robustMode = robust_mode, robustK1 = robust_k1, robustK2 = robust_k2, kSuppress = k_suppress, suppressionSigma = suppression_sigma)
    tap_fallback_used = res.converged
    base_metadata["se_tap_estimation_fallback"] = tap_fallback_used
  end
  if !res.converged
    _write_service_performance_log!(output_dir, phase_recorder, se_total_start; headline = _SE_PERF_HEADLINE, status = "failed", label = "state-estimation")
    return _api_failure("se_not_converged", "State estimation did not converge within $(max_iter) iterations.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  se_phase("postprocessing_result")

  # artifacts
  open(joinpath(output_dir, "se_diagnostics.md"), "w") do io
    # same statement the result page carries: the tap positions below are
    # model values, so the J of this run measures THEM
    tap_fallback_used && println(io, "\n> **Tap estimation fallback.** ", _SE_TAP_FALLBACK_NOTE, "\n")
    print_se_diagnostics(io, diag; topN = 15, format = :markdown)
    # topology findings of both stages append to the diagnostics report
    if !isempty(topo_pre.findings) || diag.topology_findings !== nothing
      println(io, "\n## Topology findings (advisory)\n")
      for f in topo_pre.findings
        println(io, "- precheck ", f.severity, ": `", f.kind, "` at ", f.location, " (", f.evidence, ")")
      end
      for tf in something(diag.topology_findings, NamedTuple[])
        println(io, "- classification strong: topology error suspected at station ", tf.location, " (suspects: ", tf.evidence, ")")
      end
      println(io, "\nFindings are advisory: nothing was switched or excluded. Use the hypothesis test on the result page for ranked recommendations.")
    end
  end
  view = se_view(net)
  open(joinpath(output_dir, "se_view.md"), "w") do io
    print_se_view(io, view; format = :markdown)
  end
  if res.shuntEstimates !== nothing
    open(joinpath(output_dir, "shunt_estimates.csv"), "w") do io
      println(io, "bus,bus_name,B_model_pu,B_est_pu,delta_pu,frozen")
      for r in res.shuntEstimates
        println(io, r.busIdx, ",", r.busName, ",", r.B_model, ",", r.B_est, ",", r.delta, ",", r.frozen)
      end
    end
  end
  # machine transformers are CALCULATED, never estimated (maintainer
  # directive): their tap is no state variable, so after the estimation the
  # position is back-calculated from the AVR setpoint, the dispatch P, and
  # the MEASURED machine Q (the Qinj telemetry at the machine bus when the
  # set carries it). The evaluation marks these rows as "calculated".
  tap_calc_rows = Dict{String,Any}[]
  if tap_estimation
    trm = _transformer_mrids(net)
    for k in tap_machine_idxs
      mside = _machine_side(net, k)
      mside == 0 && continue
      qrow = findfirst(mm -> mm.active && mm.typ == QinjMeas && mm.busIdx == mside, net.measurements)
      prow = findfirst(mm -> mm.active && mm.typ == PinjMeas && mm.busIdx == mside, net.measurements)
      try
        bt = calcMachineTrafoTapFromSE(net; trafo = k, p_mw = prow === nothing ? nothing : net.measurements[prow].value, q_mvar = qrow === nothing ? nothing : net.measurements[qrow].value)
        push!(tap_calc_rows, Dict{String,Any}("branch" => k, "name" => getCompName(net.branchVec[k].comp), "mrid" => get(trm, k, ""), "mode" => "ratio", "electrical_step" => round(bt.electrical_step; digits = 3), "fixed_step" => bt.fixed_step, "electrical_shift_step" => 0.0, "fixed_shift_step" => 0, "out_of_range" => false, "fixed" => false, "frozen_reason" => "none", "source" => "calculated"))
        open(logfile, "a") do io
          println(io, "machine transformer ", getCompName(net.branchVec[k].comp), ": tap CALCULATED (not estimated) from AVR setpoint and machine telemetry: electrical step ", round(bt.electrical_step; digits = 2), " -> step ", bt.fixed_step, " (Q residual ", round(bt.q_residual_mvar; digits = 2), " MVar)")
        end
      catch err
        open(logfile, "a") do io
          println(io, "machine transformer ", getCompName(net.branchVec[k].comp), ": tap back-calculation not possible (", sprint(showerror, err), ")")
        end
      end
    end
  end

  # bad data at a glance (maintainer request 2026-08-27): every suspicious
  # or eliminated measurement with its network location in one CSV; the
  # diagnostics markdown keeps the full ranking, this is the extract to
  # open first. se_state.csv stays untouched: it is the machine-read chain
  # anchor (bus voltages/balances), not a measurement report.
  suspicious = diag.diagnostics.suspicious_measurements
  # rows the robust solve acted on (final state-writing run): replacement
  # stage 3 = SUPPRESSED with the fixed sigma, staged stages 1/2 =
  # down-weighted. The bad-data extract must say so per row, not just
  # eliminated true/false.
  suppressedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage == 3)
  downweightedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage in (1, 2))
  if !isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)
    nby = _bus_name_by_idx(net)
    open(joinpath(output_dir, "se_bad_data.csv"), "w") do io
      println(io, "measurement_index,id,type,bus,from_bus,to_bus,normalized_residual,wii,localizable,eliminated,suppressed,downweighted")
      elim = Set(t.id for t in diag.eliminations)
      listed = Set{String}()
      for row in suspicious
        bus = ""
        fb = ""
        tb = ""
        if 1 <= row.measurement_index <= length(net.measurements)
          m = net.measurements[row.measurement_index]
          m.busIdx !== nothing && (bus = get(nby, m.busIdx, string(m.busIdx)))
          if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
            br = net.branchVec[m.branchIdx]
            fb = get(nby, Int(br.fromBus), string(Int(br.fromBus)))
            tb = get(nby, Int(br.toBus), string(Int(br.toBus)))
          end
        end
        push!(listed, String(row.id))
        println(io, row.measurement_index, ",", row.id, ",", row.typ, ",", bus, ",", fb, ",", tb, ",", round(row.normalized_residual; digits = 3), ",", isnan(row.wii) ? "" : round(row.wii; digits = 3), ",", row.localizable, ",", row.id in elim, ",", String(row.id) in suppressedIds, ",", String(row.id) in downweightedIds)
      end
      for t in diag.eliminations
        t.id in listed && continue
        push!(listed, String(t.id))
        println(io, t.measurement_index, ",", t.id, ",", t.typ, ",,,,", round(t.normalized_residual_before; digits = 3), ",,true,true,", String(t.id) in suppressedIds, ",", String(t.id) in downweightedIds)
      end
      # suppressed rows below the elimination threshold (possible when
      # k_suppress < k_eliminate) still belong in the extract
      for r in something(res.robustRows, NamedTuple[])
        r.stage == 3 || continue
        rid = String(r.id)
        rid in listed && continue
        push!(listed, rid)
        mi = r.measurement_index
        bus = ""
        fb = ""
        tb = ""
        ty = ""
        if 1 <= mi <= length(net.measurements)
          m = net.measurements[mi]
          ty = string(m.typ)
          m.busIdx !== nothing && (bus = get(nby, m.busIdx, string(m.busIdx)))
          if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
            br = net.branchVec[m.branchIdx]
            fb = get(nby, Int(br.fromBus), string(Int(br.fromBus)))
            tb = get(nby, Int(br.toBus), string(Int(br.toBus)))
          end
        end
        println(io, mi, ",", rid, ",", ty, ",", bus, ",", fb, ",", tb, ",", round(r.t; digits = 3), ",,,false,true,false")
      end
    end
  end
  # measured/truth/estimated delta file (generator v2): written when the
  # set carries truth_value comments (generated sets). measured and sigma
  # come from the loaded rows, truth is the generator's noise-free value,
  # estimated is measured - residual at the final diagnostics state;
  # eliminated rows keep measured/truth but no estimate contribution flag.
  # Tap rows compare the set's documented generation deviation with the
  # estimation's fixed step (authoritative per-tap table: se_tap_estimates.csv).
  deltas_written = false
  truthById = Dict{String,Float64}()
  for (k, v) in get(case_provenance, "truth_values", Dict{String,Any}())
    v isa Real && (truthById[String(k)] = Float64(v))
  end
  for line in (from_case_file ? String[] : eachline(measurement_file))
    startswith(line, "# truth_value,") || continue
    parts = split(chopprefix(line, "# truth_value,"), ",")
    length(parts) >= 2 || continue
    tv = tryparse(Float64, String(parts[end]))
    tv === nothing && continue
    truthById[String(join(parts[1:(end-1)], ","))] = tv
  end
  if !isempty(truthById)
    nbyD = _bus_name_by_idx(net)
    elimD = Set(t.id for t in diag.eliminations)
    residById = Dict{String,Float64}()
    for row in diag.final_diagnostics.measurement_ranking
      residById[String(row.id)] = row.residual
    end
    open(joinpath(output_dir, "se_deltas.csv"), "w") do io
      println(io, "# sparlectra-se-deltas v1")
      println(io, "kind,id,type,bus,from_bus,to_bus,measured,truth,estimated,delta_meas_truth,delta_est_truth,delta_est_meas,sigma,t_est,eliminated")
      fmt(v) = repr(round(v; sigdigits = 8))
      for m in net.measurements
        haskey(truthById, m.id) || continue
        tv = truthById[m.id]
        busD = m.busIdx !== nothing ? get(nbyD, m.busIdx, string(m.busIdx)) : ""
        fbD = ""
        tbD = ""
        if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
          brD = net.branchVec[m.branchIdx]
          fbD = get(nbyD, Int(brD.fromBus), string(Int(brD.fromBus)))
          tbD = get(nbyD, Int(brD.toBus), string(Int(brD.toBus)))
        end
        rres = get(residById, m.id, nothing)
        estv = rres === nothing ? nothing : m.value - rres
        println(io, "meas,", m.id, ",", m.typ, ",", busD, ",", fbD, ",", tbD, ",", repr(m.value), ",", repr(tv), ",", estv === nothing ? "" : fmt(estv), ",", fmt(m.value - tv), ",", estv === nothing ? "" : fmt(estv - tv), ",", estv === nothing ? "" : fmt(estv - m.value), ",", repr(m.sigma), ",", rres === nothing ? "" : fmt(abs(rres) / m.sigma), ",", m.id in elimD)
      end
      devByBranch = Dict{Int,Float64}(d.branch => d.steps for d in set_tap_devs)
      for t in something(res.tapEstimates, NamedTuple[])
        dev = get(devByBranch, t.branch, 0.0)
        # a tap is never measured: `measured` and the deltas against it stay
        # empty rather than repeating the truth and pretending a zero residual
        println(io, "tap,", t.name, ",", t.mode, ",,,,,", dev, ",", t.fixed_step_1, ",,", fmt(Float64(t.fixed_step_1) - dev), ",,,,false")
      end
      deltas_written = true
    end
  end

  if res.tapEstimates !== nothing || !isempty(tap_calc_rows)
    open(joinpath(output_dir, "se_tap_estimates.csv"), "w") do io
      println(io, "branch,name,mrid,mode,alpha_deg,electrical_step,fixed_step,electrical_shift_step,fixed_shift_step,r1_est,r2_est,out_of_range,fixed,frozen_reason,source")
      for t in something(res.tapEstimates, NamedTuple[])
        println(io, t.branch, ",", t.name, ",", t.mrid, ",", t.mode, ",", t.alpha_deg, ",", round(t.electrical_step_1; digits = 4), ",", t.fixed_step_1, ",", round(t.electrical_step_2; digits = 4), ",", t.fixed_step_2, ",", t.r1_est, ",", t.r2_est, ",", t.out_of_range, ",", t.fixed, ",", t.frozen_reason == :none ? "" : t.frozen_reason, ",estimated")
      end
      for c in tap_calc_rows
        println(io, c["branch"], ",", c["name"], ",", c["mrid"], ",ratio,0.0,", c["electrical_step"], ",", c["fixed_step"], ",0.0,0,,,false,false,,calculated")
      end
    end
  end
  se_phase("writing_artifacts")
  writeSEStateCSV(net; file = joinpath(output_dir, "se_state.csv"))

  # band verdict from the POST-elimination report: it must describe the
  # same state as the headline J (the pre-elimination report still tells
  # its story in se_diagnostics.md)
  verdict = diag.final_diagnostics.objective
  open(logfile, "a") do io
    println(io, "State estimation on ", basename(case_path))
    println(io, "measurements: ", summary.total, " rows from ", from_case_file ? "the case file" : basename(String(measurement_file)))
    println(io, "observability: ", obs.quality, :structural_islands in obs.notes ? " (structural islands)" : "")
    println(io, "converged in ", res.iterations, " iteration(s); J = ", round(res.objectiveJ; digits = 6), ", dof = ", res.dof, ", band reason = ", verdict.reason, res.activeObjective === nothing ? "" : string("; J_active = ", round(res.activeObjective.j; digits = 6), " (dof ", res.activeObjective.dof, ", without ", res.activeObjective.suppressed, " suppressed row(s); band verdict stays on J)"))
    if diag.topology_findings !== nothing
      # the stage-2 classification REPLACES the eliminations-exhausted
      # interpretation: this is a suspected topology error, not bad data
      for tf in diag.topology_findings
        println(io, "topology error suspected at station ", tf.location, " (eliminations exhausted, band high; suspects: ", tf.evidence, isempty(tf.notes) ? "" : string("; notes: ", join(tf.notes, ", ")), ")")
      end
      println(io, "suspicious measurements: ", length(diag.diagnostics.suspicious_measurements), "; eliminations: ", length(diag.eliminations), " (stop: ", diag.stop_reason, ", classified as topology, use the hypothesis test)")
    else
      println(io, "suspicious measurements: ", length(diag.diagnostics.suspicious_measurements), "; eliminations: ", length(diag.eliminations), " (stop: ", diag.stop_reason, ")")
    end
    effRobustMode = robust_mode === :off && robust ? :staged : robust_mode
    if effRobustMode === :staged
      println(io, "robust R modification active (staged, k1=", robust_k1, ", k2=", robust_k2, "); ", res.robustRows === nothing ? 0 : length(res.robustRows), " row(s) weighted down in the final iteration")
    elseif effRobustMode === :replacement
      println(io, "bad-data suppression active (replacement, k_suppress=", k_suppress, ", suppression sigma=", suppression_sigma, "); ", res.robustRows === nothing ? 0 : length(res.robustRows), " row(s) suppressed on the converged state; statistics keep the original sigmas")
    end
    if res.tapEstimates !== nothing
      if res.tapFixation !== nothing
        tf = res.tapFixation
        println(io, "tap estimation: ", length(res.tapEstimates), " transformer(s); J before fixation = ", round(tf.j_before; sigdigits = 4), " (dof ", tf.dof_before, "), J after fixation = ", round(tf.j_after; sigdigits = 4), " (dof ", tf.dof_after, ")", tf.fixed ? "" : "; NOT fixed (estimation run did not converge)")
        tf.offgrid_residual && println(io, "note :offgrid_tap_residual: the band failure appears only through the fixation; this is an off-grid true position or a wrong step table, NOT bad data")
      else
        println(io, "tap estimation: every released tap was frozen by the guards; no tap state entered the solve")
      end
      for t in res.tapEstimates
        println(io, "  ", t.name, isempty(t.mrid) ? "" : " (mRID " * t.mrid * ")", ": electrical step ", round(t.electrical_step_1; digits = 2), " -> fixed step ", t.fixed_step_1, t.mode == :ratio ? "" : string(", shift step ", round(t.electrical_step_2; digits = 2), " -> ", t.fixed_step_2), t.out_of_range ? " (OUT OF RANGE, clamped)" : "", t.frozen_reason == :none ? "" : " (FROZEN: $(t.frozen_reason))")
      end
    end
    (!isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)) && println(io, "bad data: se_bad_data.csv (", length(suspicious), " suspicious, ", length(diag.eliminations), " eliminated, ", length(suppressedIds), " suppressed, ", length(downweightedIds), " down-weighted, with network locations)")
    deltas_written && println(io, "deltas: se_deltas.csv (measured vs generator truth vs estimated, per row and per released tap)")
    println(io, "Artifacts: measurements.csv, se_diagnostics.md, se_view.md, se_state.csv", res.shuntEstimates !== nothing ? ", shunt_estimates.csv" : "", res.tapEstimates !== nothing ? ", se_tap_estimates.csv" : "", (!isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)) ? ", se_bad_data.csv" : "", deltas_written ? ", se_deltas.csv" : "")
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => String(format),
      "se_measurement_rows" => summary.total,
      "se_converged" => res.converged,
      # the maximum over the solves of the run, see _se_reported_iterations
      "se_iterations" => _se_reported_iterations(res),
      "se_iterations_last_solve" => res.iterations,
      "se_objective" => res.objectiveJ,
      "se_dof" => res.dof,
      "se_objective_active" => res.activeObjective === nothing ? nothing : res.activeObjective.j,
      "se_dof_active" => res.activeObjective === nothing ? nothing : res.activeObjective.dof,
      "se_suppressed_rows" => res.activeObjective === nothing ? 0 : res.activeObjective.suppressed,
      "se_band_reason" => String(verdict.reason),
      "se_j_within_3sigma" => res.jWithin3Sigma,
      "se_suspicious" => length(diag.diagnostics.suspicious_measurements),
      "se_eliminations" => length(diag.eliminations),
      "se_robust" => robust,
      "se_deltas" => deltas_written,
      "se_robust_mode" => String(robust_mode === :off && robust ? :staged : robust_mode),
      # the numbers that DECIDED this run, so a result can be read without
      # guessing which of configuration, service or form won (the reason
      # this is here at all: a run estimated with k_suppress 6.0 while the
      # configuration said 4.0, and nothing in the result showed it)
      "se_max_iter" => max_iter,
      "se_tol" => tol,
      "se_k_eliminate" => k_eliminate,
      "se_robust_k1" => robust_k1,
      "se_robust_k2" => robust_k2,
      "se_k_suppress" => k_suppress,
      "se_suppression_sigma" => suppression_sigma,
      "se_update_shunts" => update_shunts,
      "se_islands_estimated" => res.islands === nothing ? 1 : count(i -> i.estimated, res.islands),
      "se_islands_skipped" => res.islands === nothing ? 0 : count(i -> !i.estimated, res.islands),
      "se_tap_estimation" => tap_estimation,
      "se_tap_count" => res.tapEstimates === nothing ? 0 : length(res.tapEstimates),
      "se_tap_skipped_machine" => tap_skipped_machine,
      "se_tap_frozen" => res.tapEstimates === nothing ? 0 : count(t -> t.frozen_reason != :none, res.tapEstimates),
      "se_tap_fixed" => res.tapFixation === nothing ? false : res.tapFixation.fixed,
      "se_tap_j_before" => res.tapFixation === nothing ? nothing : res.tapFixation.j_before,
      "se_tap_dof_before" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_before,
      "se_tap_j_after" => res.tapFixation === nothing ? nothing : res.tapFixation.j_after,
      "se_tap_dof_after" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_after,
      "se_tap_offgrid_residual" => res.tapFixation === nothing ? false : res.tapFixation.offgrid_residual,
      "se_topology_station_findings" => diag.topology_findings === nothing ? nothing : [Dict{String,Any}("location" => tf.location, "evidence" => tf.evidence, "notes" => [String(n) for n in tf.notes]) for tf in diag.topology_findings],
      # per-trafo rows for the result page table (steps, not raw r values)
      "se_tap_estimates" => res.tapEstimates === nothing && isempty(tap_calc_rows) ? nothing : vcat([Dict{String,Any}("branch" => t.branch, "name" => t.name, "mrid" => t.mrid, "mode" => String(t.mode), "electrical_step" => round(t.electrical_step_1; digits = 3), "fixed_step" => t.fixed_step_1, "electrical_shift_step" => round(t.electrical_step_2; digits = 3), "fixed_shift_step" => t.fixed_step_2, "out_of_range" => t.out_of_range, "fixed" => t.fixed, "frozen_reason" => String(t.frozen_reason), "source" => "estimated") for t in something(res.tapEstimates, NamedTuple[])], tap_calc_rows),
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  # J/dof FIRST: J alone says nothing (it grows with the number of
  # measurements), and a reader who sees "J = 104" on a 14-bus case reads an
  # alarm where J/dof = 1.1 says the set is healthy.
  message = string(
    "State estimation completed - ", _se_reported_iterations(res), " iteration(s), J/dof = ",
    res.dof > 0 ? string(round(res.objectiveJ / res.dof; digits = 2)) : "n/a",
    " (J = ", round(res.objectiveJ; digits = 3), ", dof ", res.dof, ", band ", verdict.reason, "), ",
    length(diag.diagnostics.suspicious_measurements), " suspicious, ",
    length(diag.eliminations), " eliminated.",
    diag.topology_findings === nothing ? "" : string(" Topology error suspected at ", join((tf.location for tf in diag.topology_findings), ", "), "; run the hypothesis test."),
    res.tapFixation === nothing ? "" : string(" Tap fixation (", length(res.tapEstimates), " transformer(s)): J ", round(res.tapFixation.j_before; sigdigits = 3), " -> ", round(res.tapFixation.j_after; sigdigits = 3), res.tapFixation.offgrid_residual ? "; off-grid tap residual, not bad data" : "", "."),
    get(base_metadata, "se_set_tap_deviation", false) == true ? (isempty(get(base_metadata, "se_auto_released_taps", Int[])) ? " Hint: the set documents a tap deviation that could not be released automatically; enable 'estimate taps'." : string(" Tap estimation was released automatically on ", length(get(base_metadata, "se_auto_released_taps", Int[])), " transformer(s) named by the set.")) : "",
    duplicate_rows > 0 ? string(" Warning: ", duplicate_rows, " measurement(s) repeat an already measured quantity - that alone inflates J; regenerate the set.") : "",
  )
  # timing file plus the raw phase sequence in result.json, same as the
  # power-flow path, so the Web UI run page can show where the time went
  metadata["service_phase_timings"] = _write_service_performance_log!(output_dir, phase_recorder, se_total_start; headline = _SE_PERF_HEADLINE, label = "state-estimation")
  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = nothing,
    message = message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end

"""
    _run_pf_from_se_service(case_path, config_file, output_dir, run_id, se_state_file, se_run_id, se_mode) -> SparlectraApiResult

SE-started power-flow chain run: builds the net,
restores the estimated state from `se_state_file` (the `se_state.csv` of a
preceding SE run), and runs `runpf_from_se!` with `se_mode` (`"se_state"` or
`"se_snapshot"`). The persistent model is never mutated (the snapshot's
balance takeover lives on the working import only). Metadata records the
start mode, the referenced SE run id, and the slack pickup.
"""
function _run_pf_from_se_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, se_state_file::AbstractString, se_run_id::AbstractString, se_mode::AbstractString)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "powerflow_se_start", "se_run_id" => String(se_run_id), "se_start_mode" => String(se_mode))

  se_mode in ("se_state", "se_snapshot") || return _api_failure("invalid_request", "se_start_mode must be \"se_state\" or \"se_snapshot\", got \"$(se_mode)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)

  # the same precedence the power-flow path uses (resolve_config, D5):
  # case configuration file, the case file's deprecated block, general
  # file, defaults
  config = try
    resolve_config(config_file, case_path).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  imported, format, failure = _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata)
  failure !== nothing && return failure
  net = imported.net
  config = imported.config

  isfile(se_state_file) || return _api_failure("se_state_missing", "SE state artifact not found for run $(se_run_id); run the state estimation first.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  r = try
    readSEStateCSV!(net; file = se_state_file)
    runpf_from_se!(net, 40, 1e-8, 0; mode = Symbol(se_mode), method = :rectangular)
  catch err
    err isa PowerFlowAborted && rethrow()
    return _api_failure("se_start_failed", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  r.converged || return _api_failure("powerflow_not_converged", "SE-started power flow did not converge ($(r.iterations) iterations).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  calcNetLosses!(net)

  # deviation of the PF solution from the estimated (chosen) SE state: the
  # per-bus answer to "how far does the model pull away from the measured
  # operating point". Full table as artifact, extremes in metadata/log.
  st = _se_start_state(net)
  devrows = NamedTuple[]
  if st !== nothing
    iso = Set(net.isoNodes)
    name_by_idx = _bus_name_by_idx(net)
    for (i, nd) in enumerate(net.nodeVec)
      i in iso && continue
      vm_pf = something(nd._vm_pu, NaN)
      va_pf = something(nd._va_deg, NaN)
      dva = va_pf - st.va[i]
      # angle wrap to (-180, 180]
      dva = mod(dva + 180.0, 360.0) - 180.0
      push!(devrows, (bus = i, name = get(name_by_idx, i, string(i)), mrid = _bus_mrid(net, i; name_by_idx = name_by_idx), vm_se = st.vm[i], vm_pf = vm_pf, dvm = vm_pf - st.vm[i], va_se = st.va[i], va_pf = va_pf, dva = dva))
    end
    open(joinpath(output_dir, "se_pf_deviation.csv"), "w") do io
      println(io, "bus,name,mrid,vm_se_pu,vm_pf_pu,dvm_pu,va_se_deg,va_pf_deg,dva_deg")
      for rw in devrows
        println(io, rw.bus, ",", rw.name, ",", rw.mrid, ",", rw.vm_se, ",", rw.vm_pf, ",", rw.dvm, ",", rw.va_se, ",", rw.va_pf, ",", rw.dva)
      end
    end
  end
  maxdvm = isempty(devrows) ? NaN : maximum(abs(rw.dvm) for rw in devrows)
  maxdva = isempty(devrows) ? NaN : maximum(abs(rw.dva) for rw in devrows)
  busdvm = isempty(devrows) ? "" : devrows[argmax([abs(rw.dvm) for rw in devrows])].name
  busdva = isempty(devrows) ? "" : devrows[argmax([abs(rw.dva) for rw in devrows])].name
  meandvm = isempty(devrows) ? NaN : sum(abs(rw.dvm) for rw in devrows) / length(devrows)

  open(logfile, "a") do io
    println(io, "SE-started power flow (", se_mode, ") on ", basename(case_path))
    println(io, "start state: se_state.csv of SE run ", se_run_id)
    println(io, "converged in ", r.iterations, " iteration(s)")
    println(io, "slack pickup: ", round(r.slack_pickup_mw; digits = 6), " MW / ", round(r.slack_pickup_mvar; digits = 6), " MVar")
    if !isempty(devrows)
      println(io, "deviation PF vs SE state: max |dVm| = ", round(maxdvm; sigdigits = 4), " pu at ", busdvm, ", mean |dVm| = ", round(meandvm; sigdigits = 4), " pu, max |dVa| = ", round(maxdva; sigdigits = 4), " deg at ", busdva, " (full table: se_pf_deviation.csv)")
    end
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "iterations" => r.iterations,
      "slack_pickup_mw" => r.slack_pickup_mw,
      "slack_pickup_mvar" => r.slack_pickup_mvar,
      "se_pf_max_dvm_pu" => maxdvm,
      "se_pf_max_dvm_bus" => busdvm,
      "se_pf_mean_dvm_pu" => meandvm,
      "se_pf_max_dva_deg" => maxdva,
      "se_pf_max_dva_bus" => busdva,
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = true,
    iterations = r.iterations,
    reason = nothing,
    message = string("SE-started power flow (", se_mode, ") converged in ", r.iterations, " iteration(s); slack pickup ", round(r.slack_pickup_mw; digits = 4), " MW."),
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end
