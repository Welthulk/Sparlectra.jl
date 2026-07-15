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

function _runner_verbose(cfg::SparlectraConfig)::Int
  return cfg.output.console_diagnostics === :full ? 1 : 0
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
        controllers = collect_outer_controllers(net)
        if isempty(controllers)
          ite, status = runpf!(net, pf_cfg; verbose = verbose, pv_table_rows = qlimit_preview_rows, performance_profile = performance_profile)
          (ite, status, :none)
        else
          control_result = run_control!(net; controllers = controllers, pf_config = pf_cfg, control_config = cfg.control, verbose = verbose, performance_profile = performance_profile)
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
