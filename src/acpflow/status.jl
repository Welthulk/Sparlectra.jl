"""
    SparlectraRunResult

Typed result returned by [`run_sparlectra`](@ref). The solved network is always
available as `net`, while the remaining fields provide stable framework-level
solver, validation, control-loop, and timing metadata.
"""
struct SparlectraRunResult
  net::Net
  outcome::Symbol
  numerical_converged::Bool
  solution_available::Bool
  limit_validation_status::Symbol
  final_converged::Bool
  reason::Symbol
  reason_text::String
  iterations::Int
  elapsed_s::Float64
  solver_elapsed_s::Union{Nothing,Float64}
  final_mismatch::Float64
  method::Symbol
  control_status::Symbol
  performance_profile::Any
end

function _rectangular_solution_available(rect_status)::Bool
  rect_status === nothing && return false
  Bool(rect_status.numerical_converged) || return false

  outcome = Symbol(rect_status.status)
  reason = Symbol(rect_status.reason)
  outcome in (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented) && return false
  reason in (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented, :rescue_requested_but_not_available) && return false
  return true
end

function _rectangular_run_status(rect_status)
  status = rect_status === nothing ? :solver_error : Symbol(rect_status.status)
  outcome = status
  numerical_converged = rect_status !== nothing && Bool(rect_status.numerical_converged)
  return (
    outcome = outcome,
    numerical_converged = numerical_converged,
    solution_available = _rectangular_solution_available(rect_status),
    limit_validation_status = rect_status === nothing || !numerical_converged ? :skip : Bool(rect_status.q_limit_active_set_ok) ? :ok : :fail,
    final_converged = rect_status !== nothing && Bool(rect_status.final_converged),
    reason = rect_status === nothing ? :solver_error : Symbol(rect_status.reason),
    reason_text = rect_status === nothing ? "solver status unavailable" : String(rect_status.reason_text),
    final_mismatch = rect_status === nothing ? NaN : Float64(rect_status.final_mismatch),
  )
end

_is_accepted_control_status(status::Symbol)::Bool = status in (:none, :converged, :disabled, :no_active_controllers, :no_controllers)

function _compose_framework_status(base_status, control_status::Symbol)
  _is_accepted_control_status(control_status) && return base_status
  if control_status === :pf_failed
    meaningful_rejection = base_status.outcome !== :converged || base_status.reason !== :none || !base_status.final_converged
    meaningful_rejection && base_status.numerical_converged && return merge(base_status, (final_converged = false,))
    return merge(base_status, (outcome = :pf_failed, final_converged = false, reason = :pf_failed, reason_text = "power flow failed during controlled framework run"))
  end
  if !base_status.numerical_converged
    return merge(base_status, (outcome = :pf_failed, final_converged = false, reason = :pf_failed, reason_text = "power flow failed during controlled framework run"))
  end
  outcome = control_status === :blocked ? :control_blocked : control_status === :max_outer_iterations ? :control_max_outer_iterations : :control_not_converged
  return merge(base_status, (outcome = outcome, final_converged = false, reason = outcome, reason_text = "control loop ended with status $(control_status)"))
end

function _build_sparlectra_result(net::Net, cfg::SparlectraConfig, execution, performance_profile)::SparlectraRunResult
  status = if cfg.powerflow.method === :rectangular
    _rectangular_run_status(rectangular_pf_status(net))
  else
    converged = execution.erg == 0
    (outcome = converged ? :converged : :not_converged, numerical_converged = converged, solution_available = converged, limit_validation_status = :skip, final_converged = converged, reason = converged ? :none : :nr_mismatch_not_converged, reason_text = converged ? "none" : "NR mismatch did not converge", final_mismatch = NaN)
  end
  status = _compose_framework_status(status, execution.control_status)
  return SparlectraRunResult(net, status.outcome, status.numerical_converged, status.solution_available, status.limit_validation_status, status.final_converged, status.reason, status.reason_text, execution.iterations, execution.elapsed_s, execution.solver_elapsed_s, status.final_mismatch, cfg.powerflow.method, execution.control_status, performance_profile)
end

function _sparlectra_status_row(result::SparlectraRunResult)
  return (converged = result.final_converged, erg = result.final_converged ? 0 : 1, numerical_converged = result.numerical_converged, iterations = result.iterations, elapsed_s = result.elapsed_s, solver_elapsed_s = result.solver_elapsed_s, method = result.method, outcome = result.outcome, control_status = result.control_status, numerical_solution = result.numerical_converged ? "OK" : "FAIL", solution_available = result.solution_available, limit_validation_status = result.limit_validation_status, final_converged = result.final_converged, reason = result.reason, reason_text = result.reason_text, final_mismatch = result.final_mismatch)
end
