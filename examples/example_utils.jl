# Copyright 2023–2026 Udo Schmitz

module ExampleUtils

using Sparlectra
using LinearAlgebra

export load_yaml_config, yaml_path_from_inputs, as_runtime_thread_count,
  configure_runtime_threads!, print_runtime_thread_config

function load_yaml_config(path::AbstractString)
  isempty(path) && return Dict{String,Any}()
  isfile(path) || error("YAML config file not found: $path")
  try
    return Sparlectra.load_yaml_dict(path)
  catch
  end
  return _load_yaml_legacy(path)
end

function _parse_yaml_scalar(raw::AbstractString)
  s = strip(raw)
  isempty(s) && return nothing
  if (startswith(s, "\"") && endswith(s, "\"")) || (startswith(s, "'") && endswith(s, "'"))
    return s[2:end-1]
  end
  ls = lowercase(s)
  ls == "true" && return true
  ls == "false" && return false
  ls == "null" && return nothing
  iv = tryparse(Int, s)
  !isnothing(iv) && return iv
  fv = tryparse(Float64, s)
  !isnothing(fv) && return fv
  return s
end

function _load_yaml_legacy(path::AbstractString)
  cfg = Dict{String,Any}()
  current_section = nothing
  pending_list_key = nothing
  pending_list_section = nothing
  for rawline in eachline(path)
    uncommented = rstrip(_strip_yaml_inline_comment(rawline))
    stripped = strip(uncommented)
    isempty(stripped) && continue
    startswith(stripped, "#") && continue

    indent = firstindex(uncommented)
    while indent <= lastindex(uncommented) && uncommented[indent] == ' '
      indent = nextind(uncommented, indent)
    end
    is_nested = indent > firstindex(uncommented)

    if startswith(stripped, "-") && !isnothing(pending_list_key)
      item_raw = strip(stripped[2:end])
      target = isnothing(pending_list_section) ? cfg : cfg[pending_list_section]
      target[pending_list_key] isa AbstractVector || (target[pending_list_key] = Any[])
      push!(target[pending_list_key], _parse_yaml_scalar(item_raw))
      continue
    end

    occursin(":", stripped) || continue
    key, value = split(stripped, ":"; limit = 2)
    key = strip(key)
    value = strip(value)
    target = cfg
    if is_nested && !isnothing(current_section) && cfg[current_section] isa Dict{String,Any}
      target = cfg[current_section]
    elseif !is_nested
      current_section = nothing
    end

    if isempty(value)
      if is_nested
        target[key] = Any[]
        pending_list_key = key
        pending_list_section = current_section
      else
        cfg[key] = Dict{String,Any}()
        current_section = key
        pending_list_key = key
        pending_list_section = nothing
      end
    elseif startswith(value, "[") && endswith(value, "]")
      inner = strip(value[2:end-1])
      target[key] = isempty(inner) ? Any[] : [_parse_yaml_scalar(part) for part in split(inner, ",")]
      pending_list_key = nothing
      pending_list_section = nothing
    else
      target[key] = _parse_yaml_scalar(value)
      pending_list_key = nothing
      pending_list_section = nothing
    end
  end
  return cfg
end

function _strip_yaml_inline_comment(raw::AbstractString)
  buf = IOBuffer()
  qchar = '\0'
  for ch in raw
    if qchar != '\0'
      ch == qchar && (qchar = '\0')
      print(buf, ch)
    elseif ch == '\'' || ch == '"'
      qchar = ch
      print(buf, ch)
    elseif ch == '#'
      break
    else
      print(buf, ch)
    end
  end
  return String(take!(buf))
end

function yaml_path_from_inputs(; env_var::AbstractString = "SPARLECTRA_CONFIGURATION_YAML", fallback_paths::AbstractVector{<:AbstractString} = String[])
  !isempty(ARGS) && return ARGS[1]
  env_path = get(ENV, env_var, "")
  !isempty(env_path) && return env_path
  for p in fallback_paths
    isfile(p) && return p
  end
  return ""
end

function as_runtime_thread_count(v)
  v === nothing && return nothing
  v === false && return nothing
  v isa Integer && return Int(v)

  s = lowercase(strip(String(v)))
  s in ("", "keep", "default", "none", "false", "off", "0") && return nothing
  s == "auto" && return Sys.CPU_THREADS

  iv = tryparse(Int, s)
  if isnothing(iv) || iv < 1
    @warn "Invalid runtime thread count; keeping current setting" value = v
    return nothing
  end
  return iv
end

function configure_runtime_threads!(; julia_threads = nothing, blas_threads = nothing)
  requested_julia_threads = as_runtime_thread_count(julia_threads)
  requested_blas_threads = as_runtime_thread_count(blas_threads)

  actual_julia_threads = Threads.nthreads()
  blas_threads_before = BLAS.get_num_threads()

  if !isnothing(requested_blas_threads)
    BLAS.set_num_threads(requested_blas_threads)
  end

  blas_threads_after = BLAS.get_num_threads()
  return (;
    cpu_threads = Sys.CPU_THREADS,
    requested_julia_threads = requested_julia_threads,
    julia_threads = actual_julia_threads,
    requested_blas_threads = requested_blas_threads,
    blas_threads_before = blas_threads_before,
    blas_threads = blas_threads_after,
    julia_threads_request_applied = isnothing(requested_julia_threads) || requested_julia_threads == actual_julia_threads,
    blas_threads_request_applied = isnothing(requested_blas_threads) || requested_blas_threads == blas_threads_after,
  )
end

function print_runtime_thread_config(io::IO, status)
  println(io, "==================== Runtime thread configuration ====================")
  println(io, "CPU threads     : ", status.cpu_threads)
  println(io, "Julia threads   : ", status.julia_threads)
  if !isnothing(status.requested_julia_threads) && !status.julia_threads_request_applied
    println(io, "Julia request   : ", status.requested_julia_threads, " (not applied; start Julia with --threads or JULIA_NUM_THREADS)")
  elseif !isnothing(status.requested_julia_threads)
    println(io, "Julia request   : ", status.requested_julia_threads, " (already active)")
  end
  println(io, "BLAS threads    : ", status.blas_threads)
  if status.blas_threads_before != status.blas_threads
    println(io, "BLAS changed    : ", status.blas_threads_before, " -> ", status.blas_threads)
  elseif !isnothing(status.requested_blas_threads)
    println(io, "BLAS request    : ", status.requested_blas_threads, " (already active)")
  end
  println(io, "======================================================================\n")
  return nothing
end

end
