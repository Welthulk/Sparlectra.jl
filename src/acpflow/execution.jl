# Copyright 2023–2026 Udo Schmitz
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

# file: src/acpflow/execution.jl
# purpose: central framework execution: dispatches the configured solver
#          (rectangular NR, DC, APSLF), applies configured controllers and
#          DC seeding, and records solver timing and island diagnostics
function _runner_verbose(cfg::SparlectraConfig)::Int
  return cfg.output.console_diagnostics === :full ? 1 : 0
end

# The controllers an APSLF run leaves inactive, named in one warning: the
# outer-loop controllers (tap changers, phase shifters, machine voltage
# control) and the prosumers with a Q(U)/P(U) controller plus the
# voltage-dependent shunts. Nothing on the net is changed.
function _warn_apslf_controllers_inactive(net::Net, controllers)
  n_vdep = count(ps -> has_qu_controller(ps) || has_pu_controller(ps), net.prosumpsVec)
  n_shunt = count(sh -> sh.status != 0 && sh.model == :VoltageDependentInjection, net.shuntVec)
  (isempty(controllers) && n_vdep == 0 && n_shunt == 0) && return nothing
  @warn "power_flow.solver=apslf runs without controllers: $(length(controllers)) outer-loop controller(s) (tap changer, phase shifter, machine voltage control) are not run, $(n_vdep) Q(U)/P(U) controller(s) and $(n_shunt) voltage-dependent shunt(s) are not applied; their static setpoints are used. Set power_flow.solver=rectangular to include them." outer_controllers = length(controllers) qu_pu_controllers = n_vdep voltage_dependent_shunts = n_shunt
  return nothing
end

function _execute_sparlectra_powerflow!(net::Net, cfg::SparlectraConfig; performance_profile = nothing)
  pf_cfg = cfg.powerflow
  verbose = _runner_verbose(cfg)
  qlimit_preview_rows = cfg.output.console_q_limit_events === :summary || cfg.output.console_q_limit_events === :off ? 0 : cfg.output.console_max_rows
  phase_callback = performance_profile isa AbstractDict ? get(performance_profile, :phase_callback, phase -> nothing) : phase -> nothing
  phase_callback("solving_powerflow")
  cfg.powerflow.islands.enabled && _write_ac_island_diagnostics!(net, cfg.powerflow, performance_profile)
  iterations = 0
  erg = 2
  control_status = :none
  elapsed_s = @elapsed begin
    try
      iterations, erg, control_status = _perf_profile_time!(performance_profile, :solver_total) do
        # Controllers declared in configuration (control.controllers, issue
        # #305) are instantiated onto the net before collection; already
        # present ones are skipped so repeated runs on one net stay
        # idempotent.
        cfg.control.enabled && !isempty(cfg.control.controllers) && applyConfiguredControllers!(net, cfg.control)
        controllers = collect_outer_controllers(net)
        # Unconditional DC-seeded start (power_flow.start_mode.dc_seed_unconditional,
        # only reachable with solver=rectangular -- rejected at config time otherwise):
        # run a full standalone DC power flow first and write its angles into net as
        # the Newton-Raphson start point, mirroring rundcpf!(seed_ac_start=true) but
        # ahead of the dispatch below, so Q-limits/wrong-branch detection/diagnostics/
        # outer control all still run unchanged afterward. The candidate-selection
        # start-mode machinery (angle_mode/try_dc_start/measure_candidates) is forced
        # off for this run only, via a local pf_cfg copy, so it cannot second-guess or
        # blend away the DC seed net already carries (project_rectangular_start is a
        # pure passthrough when start_projection=false).
        # The seed is itself a flat-magnitude start (PQ buses at 1.0 pu, PV and
        # slack at their regulated magnitudes, DC angles), so the solve reads
        # it from the net with flatstart = false: under a flat start the
        # solver would otherwise rebuild 1.0 pu / 0 degrees and drop the seed
        # angles. Since #463 the flat start no longer switches the seed off,
        # so this is the path that makes `flatstart: true` plus the seed work.
        solve_pf_cfg = pf_cfg
        if pf_cfg.start_mode.dc_seed_unconditional
          _dc_seed_rectangular_angles!(net, pf_cfg; verbose = verbose, performance_profile = performance_profile)
          solve_pf_cfg = _copy_powerflow_with(pf_cfg; start_mode = _copy_start_mode_with(pf_cfg.start_mode; start_projection = false, flatstart = false))
        end
        if pf_cfg.solver === :apslf
          # APSLF solves the static injections and never evaluates a
          # controller: the outer-loop tap/PST controllers are not run and the
          # Q(U)/P(U) controllers and voltage-dependent shunts are not applied.
          # The run goes on without them and says so (it used to refuse the
          # run, which stopped every APSLF run on a MATPOWER case whose PQ
          # generator limits the importer maps to flat P(U)/Q(U) controllers)
          _warn_apslf_controllers_inactive(net, controllers)
          ite, status = _run_apslf_powerflow!(net, pf_cfg; verbose = verbose, performance_profile = performance_profile)
          (ite, status, :none)
        elseif pf_cfg.solver === :dc
          if !isempty(controllers) || has_voltage_dependent_control(net)
            throw(ArgumentError("power_flow.solver=dc does not support active controllers: neither the outer-loop tap-changer and phase-shifting transformer controllers nor the voltage-dependent Q(U)/P(U) controllers, which act inside the Newton step. A phase-shifting transformer's current, fixed angle is still represented in the DC model; only the outer control loop that adjusts it is unsupported. Disable the controllers or set power_flow.solver=rectangular."))
          end
          ite, status = _run_dc_powerflow!(net, pf_cfg; verbose = verbose, performance_profile = performance_profile)
          (ite, status, :none)
        elseif isempty(controllers)
          # the executing SparlectraConfig owns the runtime.parallel switches;
          # forward them so a non-active config still governs the island path
          ite, status = runpf!(
            net,
            solve_pf_cfg;
            verbose = verbose,
            pv_table_rows = qlimit_preview_rows,
            performance_profile = performance_profile,
            islands_parallel_enabled = cfg.runtime.parallel.enabled,
            islands_parallel_max_tasks = parallel_max_tasks(cfg.runtime.parallel),
            islands_parallel_min_work_items = cfg.runtime.parallel.min_work_items,
          )
          (ite, status, :none)
        else
          control_result = run_control!(net; controllers = controllers, pf_config = solve_pf_cfg, control_config = cfg.control, verbose = verbose, performance_profile = performance_profile)
          (control_result.last_pf_iterations, control_result.status == :pf_failed ? 1 : 0, control_result.status)
        end
      end
    catch err
      if cfg.powerflow.islands.enabled
        existing_status = rectangular_pf_status(net)
        diagnostic_status = existing_status !== nothing ? existing_status : (;
          final_mismatch = NaN,
          iterations = iterations,
          reason = :solver_error,
          status = :failed,
          stage = :pre_nr_setup,
          exception_type = nameof(typeof(err)),
          exception_message = sprint(showerror, err),
          stacktrace_top = "",
        )
        _write_ac_island_diagnostics!(net, cfg.powerflow, performance_profile; status = diagnostic_status)
      end
      rethrow()
    end
  end
  solver_elapsed_s = _solver_elapsed_from_profile(performance_profile)
  solver_elapsed_s === nothing && (solver_elapsed_s = max(0.0, Float64(elapsed_s)))
  if cfg.powerflow.islands.enabled
    status = rectangular_pf_status(net)
    _write_ac_island_diagnostics!(net, cfg.powerflow, performance_profile; status = status, iterations = iterations)
  end
  return (iterations = iterations, erg = erg, elapsed_s = elapsed_s, solver_elapsed_s = solver_elapsed_s, control_status = control_status)
end
