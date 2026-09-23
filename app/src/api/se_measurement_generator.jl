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

# file: src/api/se_measurement_generator.jl
# purpose: measurement generator behind the Web UI demo action (SE service
#          layer): MeasurementGeneratorOptions (validated at construction)
#          and _se_generate_measurement_set, which imports a case through
#          the shared SE import paths, obtains a truth state (fresh solve or
#          an adopted run state), optionally shifts transformer taps, thins
#          flow ends and passive-node rows, makes a requested number of
#          rows critical by targeted partner removal, injects seeded noise
#          and gross errors, and writes a measurement CSV v1 set with the
#          in-file case binding. Split out of run_state_estimation_service.jl
#          (which keeps the SE run itself).

"""
    MeasurementGeneratorOptions(; kwargs...)

Options of the measurement generator `_se_generate_measurement_set`. The
fields carry the same names, defaults and meaning as the former keyword
arguments of the generator; fields without a default are required keywords.
Every combination is validated at construction, so an invalid request
throws `ArgumentError` with a user-readable text BEFORE any case is imported
(the Web UI handler shows the text as its redirect message).

Fields:
- `noise`: seeded gaussian noise on the generated values (`false` gives an
  ideal set whose J is 0 by construction).
- `gross_k`, `gross_count`: bad data of `gross_k` sigma on `gross_count`
  seed-randomly drawn telemetry rows (`gross_k = 0` for none; `gross_count`
  must be at least 1).
- `tap_steps`, `tap_count`: mechanical tap deviation of `tap_steps` steps on
  up to `tap_count` seed-randomly drawn in-service transformers (`0.0` for
  none; whole steps only, the Web UI enforces this; `tap_count` at least 1).
- `include_i`: current magnitude rows on top of the power rows.
- `sigma_u_pct`, `sigma_i_pct`, `sigma_p_pct`, `sigma_q_pct`: sigmas in
  percent of reading per quantity.
- `sigma_ia_deg`: PMU current-angle sigma in degrees; `0.0` means no angle
  rows.
- `truth_source`: `:fresh_solve` (solve the configured power flow at
  reference tightness) or `:from_run` (adopt the solved voltages of run
  `run_id` under `run_root` without re-solving; needs both and is mutually
  exclusive with a tap deviation).
- `flow_ends`: `:both` keeps both flow ends per branch, `:one_balance_aware`
  keeps one, chosen by the injection telemetry at the terminal buses.
- `passive_sigma`, `passive_as_zi`: passive-node balances as protected
  zero-injection constraints (`true`) or as balance rows at `passive_sigma`
  (MW; must be positive).
- `seed`: RNG seed; the same seed regenerates the identical set.
- `critical_count`: number of rows to make critical by targeted thinning
  (`0` = off; must be zero or positive).
"""
Base.@kwdef struct MeasurementGeneratorOptions
  noise::Bool
  gross_k::Float64
  gross_count::Int = 1
  tap_steps::Float64
  tap_count::Int = 1
  include_i::Bool
  sigma_u_pct::Float64
  sigma_i_pct::Float64
  sigma_p_pct::Float64
  sigma_q_pct::Float64
  sigma_ia_deg::Float64
  truth_source::Symbol = :fresh_solve
  run_id::String = ""
  run_root::Union{Nothing,String} = nothing
  flow_ends::Symbol = :both
  passive_sigma::Float64 = 0.05
  passive_as_zi::Bool = true
  seed::Int = 42
  critical_count::Int = 0
  # the positional inner constructor is the one @kwdef's keyword constructor
  # calls, so the validation below runs for every construction path
  function MeasurementGeneratorOptions(noise, gross_k, gross_count, tap_steps, tap_count, include_i, sigma_u_pct, sigma_i_pct, sigma_p_pct, sigma_q_pct, sigma_ia_deg, truth_source, run_id, run_root, flow_ends, passive_sigma, passive_as_zi, seed, critical_count)
    opts = new(noise, gross_k, gross_count, tap_steps, tap_count, include_i, sigma_u_pct, sigma_i_pct, sigma_p_pct, sigma_q_pct, sigma_ia_deg, truth_source, run_id, run_root, flow_ends, passive_sigma, passive_as_zi, seed, critical_count)
    _validate_generator_options(opts)
    return opts
  end
end

## The argument checks of the generator, run at construction of the options:
## every message is user-readable because the Web UI shows it verbatim.
function _validate_generator_options(opts::MeasurementGeneratorOptions)::Nothing
  opts.truth_source in (:fresh_solve, :from_run) || throw(ArgumentError("truth source must be fresh_solve or from_run"))
  opts.critical_count >= 0 || throw(ArgumentError("critical measurement count must be zero or positive"))
  opts.flow_ends in (:both, :one_balance_aware) || throw(ArgumentError("flow measurements per branch must be both or one_balance_aware"))
  opts.passive_sigma > 0.0 || throw(ArgumentError("passive-node balance sigma must be positive"))
  opts.gross_count >= 1 || throw(ArgumentError("bad data count must be at least 1"))
  opts.tap_count >= 1 || throw(ArgumentError("tap deviation transformer count must be at least 1"))
  if opts.truth_source === :from_run
    # a run state is a finished snapshot: the generator never re-solves it,
    # so a tap deviation (which needs a fresh solve of the shifted model)
    # is rejected here even if the GUI lock was bypassed
    opts.tap_steps == 0.0 || throw(ArgumentError("tap deviation requires truth state 'fresh solve'; a run state is a finished snapshot and is not re-solved"))
    isempty(strip(opts.run_id)) && throw(ArgumentError("truth state 'from run' needs a run id"))
    opts.run_root isa AbstractString || throw(ArgumentError("truth state 'from run': no run directory root available"))
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
    _se_generate_measurement_set(case_path, out_path, opts::MeasurementGeneratorOptions) -> NamedTuple

Service backend of the Web UI measurement generator (the Web UI layer never
calls a solver directly). With the settings in `opts` (see
[`MeasurementGeneratorOptions`](@ref) for every field): import the case
through the shared SE import paths, optionally shift up to `tap_count`
seed-randomly selected in-service transformers by `tap_steps` MECHANICAL
steps each on their own grids (the seed decides WHICH transformers,
non-machine declared changers are drawn first), obtain the truth state
(`truth_source = :fresh_solve` solves the configured power flow at
reference tightness; `:from_run` adopts the solved voltages of run
`run_id` under `run_root` without re-solving), optionally keep only one
balance-aware flow end per branch (`flow_ends`), rewrite passive-node
balances at `passive_sigma` or as protected zero-injection constraints
(`passive_as_zi`), thin the set by targeted partner removal until
`critical_count` rows are critical, generate a measurement CSV v1 with
percent-of-reading sigmas, optional seeded noise, optional currents and
current angles, optional gross errors (`gross_k` sigma on `gross_count`
seed-randomly selected telemetry rows; protected zero-injection rows are
never corrupted), readable name-based ids, the structured tap comment
table, and the in-file case binding.

Returns `(rows, noisy, tap_note, gross_note, island_note, truth_note,
flow_note, passive_note, critical_note)` for the caller's confirmation
message. Side effect: writes `out_path`. Throws `ArgumentError` with a
user-readable text when the case has no in-service transformer for a
requested deviation, the requested run state is missing or bound to
another case, or the power flow does not converge; the option combination
itself is already validated by the constructor of `opts`.
"""
function _se_generate_measurement_set(case_path::AbstractString, out_path::AbstractString, opts::MeasurementGeneratorOptions)
  (; noise, gross_k, gross_count, tap_steps, tap_count, include_i, sigma_u_pct, sigma_i_pct, sigma_p_pct, sigma_q_pct, sigma_ia_deg, truth_source, run_id, run_root, flow_ends, passive_sigma, passive_as_zi, seed, critical_count) = opts
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
  # parallel circuits share a component name (case118 has several such
  # pairs); an id built from the name alone collided, so the second
  # circuit's rows carried the first circuit's id, and the truth-value
  # lookup by id mixed both up. A branch whose name is not unique gets
  # its branch index appended.
  branch_name_count = Dict{String,Int}()
  for br in net.branchVec
    nm = getCompName(br.comp)
    branch_name_count[nm] = get(branch_name_count, nm, 0) + 1
  end
  branch_label(k::Int) = begin
    nm = getCompName(net.branchVec[k].comp)
    get(branch_name_count, nm, 1) > 1 ? string(nm, "#", k) : nm
  end
  row_id(m) = if m.busIdx !== nothing
    string(short(m.typ), "_", get(idnby, m.busIdx, string(m.busIdx)))
  elseif m.branchIdx !== nothing
    string(short(m.typ), "_", branch_label(Int(m.branchIdx)), "_", m.direction)
  else
    m.id
  end
  for (gi2, m) in enumerate(meas)
    newid = row_id(m)
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

  # critical measurements on request (demo of issue #394): the set is
  # thinned until at least `critical_count` rows are critical, and the
  # thinning never makes the set unobservable. Every step re-reads the
  # criticality from diag(Omega) (one selected inverse, cheap since #394)
  # and removes ONE row, chosen by a targeted partner search on the dense
  # residual covariance (see _critical_anchor_partner): the anchor row that
  # needs the fewest partner removals to become critical, and its strongest
  # partner (largest normalized |Omega_ij| / sqrt(Omega_ii Omega_jj); 1 for
  # a mutually redundant pair, so one removal makes the anchor critical).
  # Above the dense-path caps the step falls back to removing the
  # redundant telemetry row with the smallest wii (the former greedy
  # sweep, which paid a trial evaluation per candidate; here the undo
  # below replaces the trial). Ties break by id, so the result is
  # deterministic. Protected ZI rows and passive balance rows stay. Rows
  # that end critical are named in the set comments, so the demo knows
  # where a gross error would hide.
  critical_note = ""
  critical_comments = String[]
  if critical_count > 0
    passiveSetC = Set(passiveBuses)
    removed_ids = String[]
    critical_ids = String[]
    # ONE Jacobian for the whole thinning: removing rows does not change
    # the rows that stay, so every step evaluates a row subset of the same
    # matrix (a rebuild per step cost minutes on case118). The observability
    # rows follow the active-measurement order the Jacobian builder uses.
    empty!(net.measurements)
    append!(net.measurements, meas)
    jac = measurement_jacobian(net)
    empty!(net.measurements)
    # the same rank decision the estimation run makes: column-normalized
    # Jacobian and the FD-aware tolerance (rank_tol_factor * jac_eps *
    # sigma_max); the matrix default tolerance called sets observable that
    # the run then refused (Web UI run e358b49e)
    H_all = _column_normalized(jac.H)
    se_tol_cfg = state_estimation_config()
    thin_tol = isempty(H_all) ? nothing : se_tol_cfg.rank_tol_factor * se_tol_cfg.jac_eps * _sigma_max(H_all)
    row_meas = Int[r.index for r in jac.rows]   # Jacobian row -> index into meas
    keep = collect(1:length(row_meas))
    removable(k) = begin
      mi = row_meas[k]
      # Jacobian rows with index 0 are link-cluster aggregates without a
      # source row: never named, never removed, never counted
      1 <= mi <= length(meas) || return false
      m = meas[mi]
      startswith(m.id, "ZI") && return false
      (m.busIdx !== nothing && m.busIdx in passiveSetC && m.typ in (PinjMeas, QinjMeas)) && return false
      return true
    end
    row_id_of(k) = 1 <= row_meas[k] <= length(meas) ? meas[row_meas[k]].id : ""
    # the full residual covariance is m x m dense: the partner search uses
    # it up to the estimator's dense-path state cap and a row cap of its
    # own, and falls back to the smallest-wii rule above
    dense_search = size(H_all, 2) <= _SE_DENSE_LINALG_MAX_N && size(H_all, 1) <= _SE_CRITICAL_SEARCH_MAX_M
    # rows whose removal broke observability once (a partner beyond the
    # rank tolerance): restored and never tried again
    blocked = Set{Int}()
    # the evaluation of the set BEFORE the last removal, so an undo costs
    # no second selected inverse
    keep_prev = keep
    obs_prev = nothing
    last_removed = 0
    obs_c = evaluate_observability_matrix(H_all[keep, :]; tol = thin_tol)
    steps = 0
    while steps <= 2 * length(row_meas) + 2
      steps += 1
      if obs_c.quality == :not_observable
        # the omega ratio and the rank tolerance are two different tests:
        # undo the removal that failed the rank test and block that row
        last_removed == 0 && break
        push!(blocked, last_removed)
        pop!(removed_ids)
        keep = keep_prev
        obs_c = obs_prev
        last_removed = 0
      end
      crit_local = Set{Int}(obs_c.numerical_critical_measurement_indices)
      crit_real = [k for k in sort!(collect(crit_local)) if 1 <= k <= length(keep) && 1 <= row_meas[keep[k]] <= length(meas)]
      critical_ids = String[row_id_of(keep[k]) for k in crit_real]
      length(crit_real) >= critical_count && break
      wii = obs_c.criticality_wii
      isempty(wii) && break
      # rows the step may remove: telemetry, not critical, not blocked
      cand_mask = Bool[!(k in crit_local) && removable(keep[k]) && !(keep[k] in blocked) for k in eachindex(keep)]
      any(cand_mask) || break
      ids_keep = String[row_id_of(keep[k]) for k in eachindex(keep)]
      H_keep = H_all[keep, :]
      pick = 0
      if dense_search
        tol_wii = _criticality_wii_tolerance(H_keep, thin_tol)
        Ω = _residual_diagnostics(H_keep, zeros(Float64, length(keep)), ones(Float64, length(keep)); need_full_omega = true).omega
        pick = _critical_anchor_partner(Ω, wii, tol_wii, cand_mask, ids_keep; max_removals = _SE_CRITICAL_SEARCH_MAX_STEPS)
      end
      if pick == 0
        # fallback (or no anchor within the step budget): the least
        # redundant removable row itself, ties by id
        cands = [(wii[k], ids_keep[k], k) for k in eachindex(keep) if cand_mask[k]]
        sort!(cands; by = c -> (c[1], c[2]))
        pick = cands[1][3]
      end
      keep_prev = keep
      obs_prev = obs_c
      last_removed = keep[pick]
      push!(removed_ids, ids_keep[pick])
      keep = [keep[k] for k in eachindex(keep) if k != pick]
      obs_c = evaluate_observability_matrix(H_all[keep, :]; tol = thin_tol)
    end
    removed_set = Set(removed_ids)
    meas = [m for m in meas if !(m.id in removed_set)]
    if length(critical_ids) >= critical_count
      critical_note = ", $(length(critical_ids)) critical row(s) after removing $(length(removed_ids)) row(s)"
    else
      critical_note = ", critical measurements: only $(length(critical_ids)) reachable (requested $(critical_count)), $(length(removed_ids)) row(s) removed"
    end
    push!(critical_comments, string("critical_target: ", critical_count, " reached: ", length(critical_ids)))
    isempty(critical_ids) || push!(critical_comments, string("critical_rows: ", join(critical_ids, " ")))
    isempty(removed_ids) || push!(critical_comments, string("critical_removed: ", join(removed_ids, " ")))
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
  append!(gencomments, critical_comments)
  tailcomments = copy(flow_comments)
  # noise-free truth values per row (the value BEFORE noise and gross
  # error): the SE run reads these comments back and writes the
  # measured/truth/estimated delta file se_deltas.csv. Passive and ZI
  # balance rows have truth 0 by construction.
  truthById = Dict{String,Float64}()
  tmeas = generateMeasurementsFromPF(net; noise = false, stddev = stddev, relativeSigma = true, includeImag = include_i, includeIa = sigma_ia_deg > 0.0)
  for tm in tmeas
    truthById[row_id(tm)] = tm.value
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
  # the bulky per-branch and per-row blocks (flow_end choices, truth values)
  # go BEHIND the data rows: the file opens with the summary, the taps table
  # and the measurements themselves, and a preview stays readable; the
  # readers skip comment lines wherever they stand
  writeMeasurementsCSV(net; file = out_path, headerComments = tapcomments, footerComments = tailcomments, busReference = _se_case_format(case_path) === :cgmes ? :mrid : :name)
  return (rows = length(net.measurements), noisy = noise, tap_note = tap_note, gross_note = gross_note, island_note = island_note, truth_note = truth_note, flow_note = flow_note, passive_note = passive_note, critical_note = critical_note)
end

## Row cap of the targeted partner search: the search reads the full m x m
## residual covariance (dense by definition; 4000 rows are 128 MB), and the
## estimator's own full-Omega cap (20000 rows, K-matrix report) is meant
## for one report, not for one matrix per removal. Sets above this cap
## thin by the smallest-wii rule instead.
const _SE_CRITICAL_SEARCH_MAX_M = 4_000

## Step budget of the per-anchor simulation: an anchor that needs more
## partner removals than this to become critical is not a target worth
## pursuing (a fully measured bus with both flow ends costs about five).
const _SE_CRITICAL_SEARCH_MAX_STEPS = 16

## Targeted partner search of the critical thinning. `Ω` is the residual
## covariance Omega = I - H (H'H)^-1 H' of the current set with unit
## weights (the projector onto the residual space, the same matrix whose
## diagonal `wii` the criticality evaluation reads; column scaling of H
## does not change it). Removing row j updates it exactly by the rank-one
## step Omega' = Omega - Omega[:, j] Omega[j, :] / Omega_jj, so
## Omega'_ii = Omega_ii (1 - r_ij^2) with the normalized covariance
## r_ij = |Omega_ij| / sqrt(Omega_ii Omega_jj): the anchor i becomes
## critical in one removal exactly when a partner with r_ij = 1 exists
## (a mutually redundant pair), and the row with the SMALLEST wii is not
## that anchor in general (the best determined row has its redundancy
## spread thinly over many partners, measured: 346 removals for 5 critical
## rows on sp_case188 when the smallest-wii row was the anchor). So every
## candidate anchor simulates the greedy partner removal (strongest
## partner first, never a row that turned critical in the simulation) on
## the updated covariance Omega^J = Omega - Omega[:, J] Omega[J, J]^-1
## Omega[J, :] up to `max_removals` steps, and the anchor with the fewest
## steps wins (ties by id). Returns the index of that anchor's FIRST
## partner, or 0 when no anchor becomes critical within the budget. The
## caller removes one row and re-evaluates; the same anchor then wins
## again with one step fewer, so the whole plan is carried out with one
## real evaluation per removal. `cand_mask` marks the rows a step may
## remove (telemetry, not critical, not blocked); anchors come from the
## same mask because a protected row is not a demo target, while it may
## still end critical as a by-product (two rows left for two states).
function _critical_anchor_partner(Ω::AbstractMatrix{Float64}, wii::Vector{Float64}, tol_wii::Float64, cand_mask::Vector{Bool}, ids::Vector{String}; max_removals::Int)::Int
  m = size(Ω, 1)
  diag0 = Float64[max(real(Ω[k, k]), 0.0) for k = 1:m]
  best_cost = typemax(Int)
  best_id = ""
  best_pick = 0
  for i = 1:m
    (cand_mask[i] && wii[i] > tol_wii) || continue
    # simulate the greedy partner removal for anchor i
    J = Int[]
    row_i = Vector{Float64}(Ω[i, :])
    dg = copy(diag0)
    cost = max_removals + 1
    for step = 1:max_removals
      wi = max(row_i[i], 0.0)
      if wi <= tol_wii
        cost = step - 1
        break
      end
      # stop early once this anchor cannot beat the best cost so far
      step > best_cost && break
      bj = 0
      br = 0.0
      for j = 1:m
        (j == i || !cand_mask[j] || dg[j] <= tol_wii || j in J) && continue
        r = abs(row_i[j]) / sqrt(wi * dg[j])
        if r > br || (r == br && bj != 0 && ids[j] < ids[bj])
          br = r
          bj = j
        end
      end
      bj == 0 && break
      push!(J, bj)
      # updated covariance after removing J: row i and the diagonal
      ΩJJ = Matrix{Float64}(Ω[J, J])
      F = cholesky(Symmetric(ΩJJ); check = false)
      issuccess(F) || break
      ΩJ = Matrix{Float64}(Ω[J, :])
      C = F \ ΩJ
      row_i = Vector{Float64}(Ω[i, :]) .- vec(transpose(C) * Vector{Float64}(Ω[J, i]))
      dg = diag0 .- vec(sum(ΩJ .* C; dims = 1))
    end
    isempty(J) && continue
    cost <= max_removals || continue
    if cost < best_cost || (cost == best_cost && ids[i] < best_id)
      best_cost = cost
      best_id = ids[i]
      best_pick = J[1]
    end
  end
  return best_pick
end
