# validate_for002_testnetz13.jl
# Revision v4: Julia 1.12 / Revise-safe builder loading.
# v3: parser accepts SubString values and is safer under VS Code/Revise reloads.
#
# Local validation example for the converted FOR001/FOR002 test network.
# Place this file in Sparlectra/examples and keep the converted files next to it,
# or pass explicit paths via command-line flags.
#
# Defaults expected next to this file:
#   FOR002.DAT
#   case_for001_testnetz13.m
#   build_for001_testnetz13.jl
#
# Example:
#   julia --project=. examples/validate_for002_testnetz13.jl
#   julia --project=. examples/validate_for002_testnetz13.jl --contingencies
#   julia --project=. examples/validate_for002_testnetz13.jl --for002=C:\\tmp\\FOR002.DAT --matpower=C:\\tmp\\case_for001_testnetz13.m --builder=C:\\tmp\\build_for001_testnetz13.jl
# file: examples/validate_for002_testnetz13.jl
# Date: 2026-07-01

using Sparlectra
using Dates
using Printf

const DEFAULT_OUTDIR = joinpath(@__DIR__, "_out", "for002_validation")
const NUMBER_PAT = "[-+]?\\d+(?:\\.\\d*)?(?:[Ee][-+]?\\d+)?"
const NUMBER_RE = Regex(NUMBER_PAT)
const LOSS_RE = Regex("GESAMTE NETZVERLUSTE.*P\\s*=\\s*(" * NUMBER_PAT * ")\\s*MW\\s*Q\\s*=\\s*(" * NUMBER_PAT * ")\\s*MVAR")

mutable struct For002Scenario
  name::String
  outage_letter::Union{Nothing,String}
  outage_from::Union{Nothing,String}
  outage_to::Union{Nothing,String}
  buses::Dict{String,NamedTuple}
  branches::Vector{NamedTuple}
  total_p_loss_MW::Union{Nothing,Float64}
  total_q_loss_MVar::Union{Nothing,Float64}
  iterations::Union{Nothing,Int}
end

_maybe_string(x) = x === nothing ? nothing : String(x)

function For002Scenario(name::AbstractString; outage_letter = nothing, outage_from = nothing, outage_to = nothing)
  return For002Scenario(String(name), _maybe_string(outage_letter), _maybe_string(outage_from), _maybe_string(outage_to), Dict{String,NamedTuple}(), NamedTuple[], nothing, nothing, nothing)
end

struct CliOptions
  for002_path::String
  matpower_path::String
  builder_path::String
  builder_function::String
  output_dir::String
  run_matpower::Bool
  run_builder::Bool
  run_contingencies::Bool
  maxiter::Int
  solver_tol::Float64
  verbose::Int
  autodamp::Bool
  opt_flatstart::Bool
  qlimits_enabled::Bool
  tol_v_kV::Float64
  tol_vm_pu::Float64
  tol_va_deg::Float64
  tol_p_MW::Float64
  tol_q_MVar::Float64
  tol_branch_p_MW::Float64
  tol_branch_q_MVar::Float64
  branch_b_scale::Float64
  line_b_scale::Float64
  trafo_b_scale::Float64
  matpower_ratio::Symbol
  matpower_pq_gen_controllers::Bool
  builder_pq_gen_controllers::Bool
  sweep_branch_b::Bool
  infer_pv_buses::Symbol
  diagnose_active_flow::Bool
  diagnose_worst_branches::Bool
  sweep_worst_branch_angle::Bool
end

struct ModelRunResult
  model_name::String
  scenario_name::String
  branch_index_out::Union{Nothing,Int}
  net::Net
  report::ACPFlowReport
  iterations::Int
  status::Int
  residual_inf::Union{Nothing,Float64}
end

_case_tag(s::AbstractString) = replace(splitext(basename(s))[1], r"[^A-Za-z0-9_.-]" => "_")
_norm_name(s::AbstractString) = uppercase(strip(replace(String(s), r"\s+" => " ")))
function _normalize_for002_bus_name(name::AbstractString)::String
  s = String(strip(name))
  return endswith(s, " S") ? String(strip(s[1:end-2])) : s
end
_angle_delta_deg(calc_deg::Real, ref_deg::Real)::Float64 = mod(Float64(calc_deg) - Float64(ref_deg) + 180.0, 360.0) - 180.0
_csv_quote(s::AbstractString) = "\"" * replace(String(s), "\"" => "\"\"") * "\""
_csv_value(x) = x === nothing ? "" : x isa AbstractString ? _csv_quote(x) : x isa Missing ? "" : string(x)
_csv_line(xs...) = join((_csv_value(x) for x in xs), ",")

function _field(line::AbstractString, lo::Int, hi::Int)::String
  lo <= hi || return ""
  n = lastindex(line)
  lo > n && return ""
  return String(SubString(line, lo, min(hi, n)))
end

function _numbers(s::AbstractString)::Vector{Float64}
  return [parse(Float64, m.match) for m in eachmatch(NUMBER_RE, s)]
end

function _parse_cli_args(args::Vector{String})::CliOptions
  defaults = Dict{String,String}(
    "for002" => joinpath(@__DIR__, "FOR002.DAT"),
    "matpower" => joinpath(@__DIR__, "case_for001_testnetz13.m"),
    "builder" => joinpath(@__DIR__, "build_for001_testnetz13.jl"),
    "builder-function" => "build_for001_testnetz13",
    "output-dir" => DEFAULT_OUTDIR,
    "maxiter" => "60",
    "solver-tol" => "1e-8",
    "verbose" => "0",
    "tol-v-kv" => "0.25",
    "tol-vm-pu" => "0.002",
    "tol-va-deg" => "0.10",
    "tol-p-mw" => "0.75",
    "tol-q-mvar" => "1.50",
    "tol-branch-p-mw" => "1.00",
    "tol-branch-q-mvar" => "2.00",
    "branch-b-scale" => "1.0",
    "matpower-ratio" => "normal",
    "matpower-pq-gen-controllers" => "true",
    "builder-pq-gen-controllers" => "true",
    "infer-pv-buses" => "off",
  )
  line_b_scale_value = nothing
  trafo_b_scale_value = nothing
  run_matpower = true
  run_builder = true
  run_contingencies = false
  autodamp = true
  opt_flatstart = false
  qlimits_enabled = false
  sweep_branch_b = false
  diagnose_active_flow = false
  diagnose_worst_branches = false
  sweep_worst_branch_angle = false

  for arg in args
    if arg == "--contingencies"
      run_contingencies = true
    elseif arg == "--sweep-branch-b"
      sweep_branch_b = true
    elseif arg == "--diagnose-active-flow"
      diagnose_active_flow = true
    elseif arg == "--diagnose-worst-branches"
      diagnose_worst_branches = true
    elseif arg == "--sweep-worst-branch-angle"
      sweep_worst_branch_angle = true
    elseif arg == "--no-matpower"
      run_matpower = false
    elseif arg == "--no-builder"
      run_builder = false
    elseif arg == "--no-autodamp"
      autodamp = false
    elseif arg == "--flatstart"
      opt_flatstart = true
    elseif arg == "--qlimits"
      qlimits_enabled = true
    elseif arg in ("-h", "--help")
      println("Usage: julia --project=. examples/validate_for002_testnetz13.jl [options]")
      println()
      println("Options:")
      println("  --for002=PATH              FOR002.DAT reference report")
      println("  --matpower=PATH            converted MATPOWER case")
      println("  --builder=PATH             converted native Sparlectra builder file")
      println("  --builder-function=NAME    builder function name, default build_for001_testnetz13")
      println("  --output-dir=PATH          output directory")
      println("  --contingencies            also run FOR002 outage scenarios")
      println("  --no-matpower              skip MATPOWER model path")
      println("  --no-builder               skip native Sparlectra builder path")
      println("  --flatstart                use flat start for the solver")
      println("  --no-autodamp              disable autodamp")
      println("  --qlimits                  enable q-limit handling")
      println("  --maxiter=N                Newton iteration limit")
      println("  --solver-tol=X             solver tolerance")
      println("  --verbose=N                Sparlectra solver verbosity")
      println("  --tol-v-kv=X               voltage kV comparison tolerance")
      println("  --tol-vm-pu=X              voltage pu comparison tolerance")
      println("  --tol-va-deg=X             angle comparison tolerance")
      println("  --tol-p-mw=X               bus/loss active-power tolerance")
      println("  --tol-q-mvar=X             bus/loss reactive-power tolerance")
      println("  --tol-branch-p-mw=X        branch active-power tolerance")
      println("  --tol-branch-q-mvar=X      branch reactive-power tolerance")
      println("  --branch-b-scale=X         diagnostic multiplier for branch b_pu before solving")
      println("  --line-b-scale=X           diagnostic multiplier for AC-line branch b_pu before solving")
      println("  --trafo-b-scale=X          diagnostic multiplier for transformer branch b_pu before solving")
      println("  --sweep-branch-b           run compact base-case line/transformer branch-charging sweep")
      println("  --diagnose-active-flow     write active-power and branch/tap diagnostic CSVs")
      println("  --diagnose-worst-branches  write detailed top branch-flow offender diagnostics")
      println("  --sweep-worst-branch-angle diagnostic one-at-a-time phase-shift sweep for top branch-flow offenders")
      println("  --infer-pv-buses=off|from-q-limits|from-gen-vg|all-generator-buses")
      println("                              diagnostic PV-bus inference mode; default off")
      println("  --matpower-ratio=normal|reciprocal")
      println("                              MATPOWER import transformer-ratio convention for this validator")
      println("  --matpower-pq-gen-controllers=true|false")
      println("                              enable MATPOWER-import P/Q controllers on PQ generators")
      println("  --builder-pq-gen-controllers=true|false")
      println("                              keep or clear native-builder P/Q controllers on PQ generators")
      exit(0)
    elseif startswith(arg, "--") && occursin("=", arg)
      key, value = split(arg[3:end], "="; limit = 2)
      if key == "line-b-scale"
        line_b_scale_value = value
      elseif key == "trafo-b-scale"
        trafo_b_scale_value = value
      else
        haskey(defaults, key) || throw(ArgumentError("Unsupported option --$key"))
        defaults[key] = value
      end
    elseif startswith(arg, "--")
      throw(ArgumentError("Unsupported flag: $arg"))
    else
      throw(ArgumentError("Unsupported positional argument: $arg"))
    end
  end

  run_matpower || run_builder || throw(ArgumentError("At least one model path must be enabled."))
  branch_b_scale = parse(Float64, defaults["branch-b-scale"])
  line_b_scale = line_b_scale_value === nothing ? branch_b_scale : parse(Float64, line_b_scale_value)
  trafo_b_scale = trafo_b_scale_value === nothing ? branch_b_scale : parse(Float64, trafo_b_scale_value)
  matpower_ratio = Symbol(defaults["matpower-ratio"])
  matpower_ratio in (:normal, :reciprocal) || throw(ArgumentError("Unsupported --matpower-ratio=$(defaults["matpower-ratio"]); expected normal or reciprocal."))
  infer_pv_buses = Symbol(defaults["infer-pv-buses"])
  infer_pv_buses in (:off, Symbol("from-q-limits"), Symbol("from-gen-vg"), Symbol("all-generator-buses")) || throw(ArgumentError("Unsupported --infer-pv-buses=$(defaults["infer-pv-buses"]); expected off, from-q-limits, from-gen-vg, or all-generator-buses."))
  matpower_pq_gen_controllers = _parse_bool_option("--matpower-pq-gen-controllers", defaults["matpower-pq-gen-controllers"])
  builder_pq_gen_controllers = _parse_bool_option("--builder-pq-gen-controllers", defaults["builder-pq-gen-controllers"])

  return CliOptions(
    defaults["for002"],
    defaults["matpower"],
    defaults["builder"],
    defaults["builder-function"],
    defaults["output-dir"],
    run_matpower,
    run_builder,
    run_contingencies,
    parse(Int, defaults["maxiter"]),
    parse(Float64, defaults["solver-tol"]),
    parse(Int, defaults["verbose"]),
    autodamp,
    opt_flatstart,
    qlimits_enabled,
    parse(Float64, defaults["tol-v-kv"]),
    parse(Float64, defaults["tol-vm-pu"]),
    parse(Float64, defaults["tol-va-deg"]),
    parse(Float64, defaults["tol-p-mw"]),
    parse(Float64, defaults["tol-q-mvar"]),
    parse(Float64, defaults["tol-branch-p-mw"]),
    parse(Float64, defaults["tol-branch-q-mvar"]),
    branch_b_scale,
    line_b_scale,
    trafo_b_scale,
    matpower_ratio,
    matpower_pq_gen_controllers,
    builder_pq_gen_controllers,
    sweep_branch_b,
    infer_pv_buses,
    diagnose_active_flow,
    diagnose_worst_branches,
    sweep_worst_branch_angle,
  )
end

function _parse_bool_option(name::AbstractString, value::AbstractString)::Bool
  normalized = lowercase(strip(String(value)))
  normalized == "true" && return true
  normalized == "false" && return false
  throw(ArgumentError("Unsupported $name=$value; expected true or false."))
end

function _try_parse_for002_outage(line::AbstractString)
  occursin("AUSFALL DES ZWEIGES", line) || return nothing
  m = match(r"AUSFALL DES ZWEIGES\s+(?:([A-Z])\s+)?VON\s+(.+?)\s+NACH\s+(.+?)\s+I", line)
  m === nothing && return (letter = nothing, from = nothing, to = nothing, title = strip(line))
  letter = m.captures[1] === nothing ? nothing : String(strip(m.captures[1]))
  from_bus = _normalize_for002_bus_name(m.captures[2])
  to_bus = _normalize_for002_bus_name(m.captures[3])
  title = isnothing(letter) ? "outage $from_bus -> $to_bus" : "outage $letter $from_bus -> $to_bus"
  return (letter = letter, from = from_bus, to = to_bus, title = title)
end

function _try_parse_for002_bus(line::AbstractString)
  startswith(line, "I ") || return nothing
  name = _normalize_for002_bus_name(_field(line, 3, 13))
  isempty(name) && return nothing
  v = _numbers(_field(line, 15, 28))
  gen = _numbers(_field(line, 30, 44))
  load = _numbers(_field(line, 46, 60))
  length(v) == 2 && length(gen) == 2 && length(load) == 2 || return nothing
  return (bus_name = name, v_kV = v[1], va_deg = v[2], p_gen_MW = gen[1], q_gen_MVar = gen[2], p_load_MW = load[1], q_load_MVar = load[2])
end

function _try_parse_for002_branch(line::AbstractString, current_bus::Union{Nothing,AbstractString})
  current_bus === nothing && return nothing
  startswith(line, "I ") || return nothing
  # A branch row has an empty bus column and a filled branch-target field.
  strip(_field(line, 3, 13)) == "" || return nothing
  target_field = _field(line, 62, 85)
  to_bus = _normalize_for002_bus_name(_field(target_field, 6, 13))
  isempty(to_bus) && return nothing
  nr = String(strip(_field(target_field, 15, 15)))
  branch_type = String(strip(_field(target_field, 19, 22)))
  flow = _numbers(_field(line, 87, 102))
  losses = _numbers(_field(line, 104, 119))
  current = _numbers(_field(line, 121, 127))
  length(flow) == 2 || return nothing
  length(losses) == 2 || return nothing
  return (from_bus = _normalize_for002_bus_name(current_bus), to_bus = to_bus, nr = nr, type = branch_type, p_MW = flow[1], q_MVar = flow[2], p_loss_MW = losses[1], q_loss_MVar = losses[2], loading_percent = isempty(current) ? missing : current[1])
end

function parse_for002(path::AbstractString)::Vector{For002Scenario}
  isfile(path) || throw(ArgumentError("FOR002 reference file not found: $path"))
  scenarios = For002Scenario[]
  current = For002Scenario("base")
  push!(scenarios, current)
  current_bus::Union{Nothing,String} = nothing

  for raw_line in eachline(path)
    line = replace(raw_line, '\f' => ' ')

    outage = _try_parse_for002_outage(line)
    if outage !== nothing
      current = For002Scenario(outage.title; outage_letter = outage.letter, outage_from = outage.from, outage_to = outage.to)
      push!(scenarios, current)
      current_bus = nothing
      continue
    end

    bus = _try_parse_for002_bus(line)
    if bus !== nothing
      current.buses[_norm_name(bus.bus_name)] = bus
      current_bus = String(bus.bus_name)
      continue
    end

    branch = _try_parse_for002_branch(line, current_bus)
    branch !== nothing && push!(current.branches, branch)

    loss_match = match(LOSS_RE, line)
    if loss_match !== nothing
      current.total_p_loss_MW = parse(Float64, loss_match.captures[1])
      current.total_q_loss_MVar = parse(Float64, loss_match.captures[2])
      continue
    end

    iter_match = match(r"NEWTONITERATION\s+(\d+)\s+ITERATIONEN", line)
    if iter_match !== nothing
      current.iterations = parse(Int, iter_match.captures[1])
      continue
    end
  end

  return [s for s in scenarios if !isempty(s.buses) || !isempty(s.branches) || s.total_p_loss_MW !== nothing]
end

function _parse_matpower_string_vector(path::AbstractString, field_name::AbstractString)::Vector{String}
  isfile(path) || return String[]
  values = String[]
  inside = false
  start_pat = "mpc.$field_name = {"
  for line in eachline(path)
    if !inside && occursin(start_pat, line)
      inside = true
      continue
    end
    inside || continue
    occursin("};", line) && break
    m = match(r"'([^']*)'", line)
    m !== nothing && push!(values, m.captures[1])
  end
  return values
end

function _branch_nr_from_label(label::AbstractString)::String
  m = match(r"^[LT]\d+([A-Z-])\s+", String(label))
  if m === nothing
    return ""
  end
  nr = m.captures[1]
  return nr == "-" ? "" : nr
end

function _branch_key(from_bus::AbstractString, to_bus::AbstractString, nr)::Tuple{String,String,String}
  nr_s = nr === nothing || nr isa Missing ? "" : strip(String(nr))
  return (_norm_name(from_bus), _norm_name(to_bus), uppercase(nr_s))
end

function _reference_branch_dict(s::For002Scenario)::Dict{Tuple{String,String,String},NamedTuple}
  d = Dict{Tuple{String,String,String},NamedTuple}()
  for row in s.branches
    d[_branch_key(row.from_bus, row.to_bus, row.nr)] = row
  end
  return d
end

function _find_branch_index_for_outage(s::For002Scenario, branch_labels::Vector{String}, report::Union{Nothing,ACPFlowReport} = nothing)::Union{Nothing,Int}
  s.outage_from === nothing && return nothing
  want_from = _norm_name(s.outage_from)
  want_to = _norm_name(s.outage_to === nothing ? "" : s.outage_to)
  want_nr = s.outage_letter === nothing ? "" : uppercase(strip(s.outage_letter))

  # Preferred path: branch label metadata from the MATPOWER conversion.
  candidates = Int[]
  for (idx, label) in enumerate(branch_labels)
    nr = _branch_nr_from_label(label)
    has_direction = occursin("->", label)
    if has_direction
      parts = split(label, "->"; limit = 2)
      lhs = strip(parts[1])
      rhs = strip(parts[2])
      lhs_match = match(r"^[LT]\d+[A-Z-]?\s+(.+)$", lhs)
      from_label = lhs_match === nothing ? lhs : lhs_match.captures[1]
      to_label = rhs
      if _norm_name(from_label) == want_from && _norm_name(to_label) == want_to && uppercase(nr) == want_nr
        push!(candidates, idx)
      end
    end
  end
  length(candidates) == 1 && return candidates[1]

  # Fallback: use branch report direction and label-derived parallel id.
  if report !== nothing
    candidates = Int[]
    for br in report.branches
      idx = br.branch_index
      label = idx <= length(branch_labels) ? branch_labels[idx] : br.branch
      nr = _branch_nr_from_label(label)
      # report.branches carries bus indices, not names, so this fallback is intentionally limited.
      uppercase(nr) == want_nr || continue
      push!(candidates, idx)
    end
    length(candidates) == 1 && return candidates[1]
  end

  return nothing
end

function _build_v_pf_from_net(net::Net, model::PFModel)::Vector{ComplexF64}
  V = Vector{ComplexF64}(undef, length(model.busIdx_net))
  for (k, bus_idx) in enumerate(model.busIdx_net)
    node = net.nodeVec[bus_idx]
    V[k] = ComplexF64(node._vm_pu * cosd(node._va_deg), node._vm_pu * sind(node._va_deg))
  end
  return V
end

function _safe_residual_inf(net::Net)::Union{Nothing,Float64}
  try
    model = buildPfModel(net; flatstart = false, include_limits = false, verbose = 0)
    V = _build_v_pf_from_net(net, model)
    r = mismatchInf(model, V)
    return isfinite(r) ? r : nothing
  catch err
    @warn "Could not compute post-solve mismatchInf" exception = (typeof(err), sprint(showerror, err))
    return nothing
  end
end

function _branch_kind_from_metadata(net::Net, branch_index::Int)
  metadata = get(net.matpower_branch_metadata, branch_index, nothing)
  metadata === nothing && return nothing
  hasproperty(metadata, :orig_kind) || return nothing
  kind = metadata.orig_kind
  kind isa Symbol && return kind
  return Symbol(lowercase(strip(String(kind))))
end

function _branch_diagnostic_kind(net::Net, br)::Symbol
  metadata_kind = _branch_kind_from_metadata(net, br.branchIdx)
  metadata_kind in (:line, :transformer) && return metadata_kind
  br.ratio == 0.0 && return :line
  return :transformer
end

function _branch_model_kind_label(net::Net, br)::String
  kind = _branch_diagnostic_kind(net, br)
  kind === :line && return "line"
  kind === :transformer && return "transformer"
  return String(kind)
end

function _scale_branch_b!(net::Net, opt::CliOptions)
  counts = Dict(:line => 0, :transformer => 0, :unknown => 0)
  for br in net.branchVec
    kind = _branch_diagnostic_kind(net, br)
    if kind === :line
      br.b_pu *= opt.line_b_scale
      counts[:line] += 1
    elseif kind === :transformer
      br.b_pu *= opt.trafo_b_scale
      counts[:transformer] += 1
    else
      counts[:unknown] += 1
    end
  end
  return counts
end

function _scale_branch_b!(net::Net, scale::Float64)
  scale == 1.0 && return nothing
  for br in net.branchVec
    br.b_pu *= scale
  end
  return nothing
end

function _clear_pq_gen_controllers!(net::Net)
  cleared = 0
  for ps in net.prosumpsVec
    Sparlectra.isGenerator(ps) || continue
    ps.isRegulated && continue
    had_controller = Sparlectra.has_pu_controller(ps) || Sparlectra.has_qu_controller(ps)
    ps.puController = nothing
    ps.quController = nothing
    cleared += had_controller ? 1 : 0
  end
  return cleared
end

function _branch_b_diagnostic_rows(opt::CliOptions, model_builders::Vector{Tuple{String,Function}}, branch_labels::Vector{String})
  mpc = Sparlectra.MatpowerIO.read_case(opt.matpower_path; legacy_compat = true)
  bus_base_kv = Dict(Int(row[1]) => Float64(row[10]) for row in eachrow(mpc.bus))
  matpower_br_b = [Float64(row[5]) for row in eachrow(mpc.branch)]
  raw_for001_b_si = Float64[]
  for row in eachrow(mpc.branch)
    fbus = Int(row[1])
    base_kv = get(bus_base_kv, fbus, 0.0)
    br_b = Float64(row[5])
    push!(raw_for001_b_si, base_kv == 0.0 ? NaN : br_b * Float64(mpc.baseMVA) / base_kv^2)
  end

  b_by_model = Dict{String,Vector{Float64}}()
  kind_by_model = Dict{String,Vector{String}}()
  bus_by_model = Dict{String,Vector{NamedTuple}}()
  for (model_name, buildfn) in model_builders
    net = buildfn()
    _scale_branch_b!(net, opt)
    b_by_model[String(model_name)] = [br.b_pu for br in net.branchVec]
    kind_by_model[String(model_name)] = [_branch_model_kind_label(net, br) for br in net.branchVec]
    names_by_idx = _bus_name_by_idx(net)
    bus_by_model[String(model_name)] = [
      (
        from_bus = get(names_by_idx, br.fromBus, string(br.fromBus)),
        to_bus = get(names_by_idx, br.toBus, string(br.toBus)),
        from_base_kv = net.nodeVec[br.fromBus].comp.cVN,
        to_base_kv = net.nodeVec[br.toBus].comp.cVN,
      ) for br in net.branchVec
    ]
  end

  nbranch = length(matpower_br_b)
  rows = NamedTuple[]
  for idx in 1:nbranch
    matpower_net_b = get(b_by_model, "matpower", Float64[])
    builder_net_b = get(b_by_model, "builder", Float64[])
    matpower_kind = get(kind_by_model, "matpower", String[])
    builder_kind = get(kind_by_model, "builder", String[])
    bus_rows = get(bus_by_model, "matpower", get(bus_by_model, "builder", NamedTuple[]))
    bus_row = idx <= length(bus_rows) ? bus_rows[idx] : (from_bus = "", to_bus = "", from_base_kv = missing, to_base_kv = missing)
    model_kind = idx <= length(matpower_kind) ? matpower_kind[idx] : idx <= length(builder_kind) ? builder_kind[idx] : "missing"
    kind = mpc.branch_kind === nothing || idx > length(mpc.branch_kind) ? missing : mpc.branch_kind[idx]
    scale = model_kind == "line" ? opt.line_b_scale : model_kind == "transformer" ? opt.trafo_b_scale : opt.branch_b_scale
    matpower_b_scaled = matpower_br_b[idx] * scale
    sparlectra_b = idx <= length(matpower_net_b) ? matpower_net_b[idx] : idx <= length(builder_net_b) ? builder_net_b[idx] : missing
    push!(
      rows,
      (
        branch_index = idx,
        branch_label = idx <= length(branch_labels) ? branch_labels[idx] : string(idx),
        from_bus = bus_row.from_bus,
        to_bus = bus_row.to_bus,
        from_base_kv = bus_row.from_base_kv,
        to_base_kv = bus_row.to_base_kv,
        branch_kind_metadata = kind,
        sparlectra_model_kind = model_kind,
        raw_or_inferred_for001_b_S = raw_for001_b_si[idx],
        matpower_br_b_total_pu_unscaled = matpower_br_b[idx],
        matpower_br_b_total_pu_scaled = matpower_b_scaled,
        matpower_net_b_pu_after_scaling = idx <= length(matpower_net_b) ? matpower_net_b[idx] : missing,
        builder_net_b_pu_after_scaling = idx <= length(builder_net_b) ? builder_net_b[idx] : missing,
        sparlectra_per_end_b_pu = sparlectra_b === missing ? missing : sparlectra_b / 2.0,
        line_b_scale = opt.line_b_scale,
        trafo_b_scale = opt.trafo_b_scale,
        branch_b_scale = opt.branch_b_scale,
        matpower_ratio = String(opt.matpower_ratio),
        matpower_pq_gen_controllers = opt.matpower_pq_gen_controllers,
        builder_pq_gen_controllers = opt.builder_pq_gen_controllers,
      ),
    )
  end
  return rows
end

function _write_branch_b_diagnostic_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "branch_index,branch_label,from_bus,to_bus,from_base_kv,to_base_kv,branch_kind_metadata,sparlectra_model_kind,raw_or_inferred_for001_b_S,matpower_br_b_total_pu_unscaled,matpower_br_b_total_pu_scaled,matpower_net_b_pu_after_scaling,builder_net_b_pu_after_scaling,sparlectra_per_end_b_pu_after_scaling,line_b_scale,trafo_b_scale,branch_b_scale,matpower_ratio,matpower_pq_gen_controllers,builder_pq_gen_controllers")
    for row in rows
      println(
        io,
        _csv_line(
          row.branch_index,
          row.branch_label,
          row.from_bus,
          row.to_bus,
          row.from_base_kv,
          row.to_base_kv,
          row.branch_kind_metadata,
          row.sparlectra_model_kind,
          row.raw_or_inferred_for001_b_S,
          row.matpower_br_b_total_pu_unscaled,
          row.matpower_br_b_total_pu_scaled,
          row.matpower_net_b_pu_after_scaling,
          row.builder_net_b_pu_after_scaling,
          row.sparlectra_per_end_b_pu,
          row.line_b_scale,
          row.trafo_b_scale,
          row.branch_b_scale,
          row.matpower_ratio,
          row.matpower_pq_gen_controllers,
          row.builder_pq_gen_controllers,
        ),
      )
    end
  end
  return path
end

function _print_branch_b_diagnostic_table(rows::Vector{NamedTuple})
  isempty(rows) && return nothing
  println("Branch charging diagnostic table")
  println("  raw FOR001 B is inferred from converted BR_B, baseMVA, and from-side base kV.")
  println("  MATPOWER BR_B and Sparlectra Branch.b_pu are total branch charging; Sparlectra stamps b/2 at each end.")
  @printf(
    "%-4s %-34s %14s %-12s %16s %-12s %16s %16s %-12s %16s %-12s\n",
    "idx",
    "branch",
    "FOR001_B_S",
    "src_basis",
    "MATPOWER_BR_B",
    "mp_basis",
    "mp_net_b_pu",
    "builder_b_pu",
    "sp_basis",
    "per_end_b_pu",
    "stamp_basis",
  )
  for row in rows
    @printf(
      "%-4d %-34s %14.8g %-12s %16.8g %-12s %16s %16s %-12s %16s %-12s\n",
      row.branch_index,
      row.branch_label,
      row.raw_or_inferred_for001_b_S,
      "inferred",
      row.matpower_br_b_total_pu_scaled,
      "total",
      row.matpower_net_b_pu_after_scaling === missing ? "missing" : @sprintf("%.8g", row.matpower_net_b_pu_after_scaling),
      row.builder_net_b_pu_after_scaling === missing ? "missing" : @sprintf("%.8g", row.builder_net_b_pu_after_scaling),
      "total",
      row.sparlectra_per_end_b_pu === missing ? "missing" : @sprintf("%.8g", row.sparlectra_per_end_b_pu),
      "b/2/end",
    )
  end
  println()
  return nothing
end

function _run_powerflow!(net::Net, opt::CliOptions; branch_index_out::Union{Nothing,Int} = nothing)
  if branch_index_out !== nothing
    setNetBranchStatus!(net = net, branchNr = branch_index_out, status = 0)
  end
  iters, status = runpf!(net, opt.maxiter, opt.solver_tol, opt.verbose; method = :rectangular, autodamp = opt.autodamp, opt_flatstart = opt.opt_flatstart, qlimits_enabled = opt.qlimits_enabled, wrong_branch_detection = :warn)
  calcNetLosses!(net)
  residual = _safe_residual_inf(net)
  report = buildACPFlowReport(net; ite = iters, tol = opt.solver_tol, converged = (status == 0), solver = :NR)
  return report, iters, status, residual
end

function _run_model(model_name::AbstractString, buildfn::Function, opt::CliOptions, scenario_name::AbstractString; branch_index_out::Union{Nothing,Int} = nothing)::ModelRunResult
  net = buildfn()
  _scale_branch_b!(net, opt)
  report, iters, status, residual = _run_powerflow!(net, opt; branch_index_out = branch_index_out)
  return ModelRunResult(String(model_name), String(scenario_name), branch_index_out, net, report, iters, status, residual)
end

function _bus_result_dict(report::ACPFlowReport)::Dict{String,NamedTuple}
  d = Dict{String,NamedTuple}()
  for row in report.nodes
    d[_norm_name(row.bus_name)] = row
  end
  return d
end

function _bus_name_by_idx(net::Net)::Dict{Int,String}
  d = Dict{Int,String}()
  for (name, idx) in net.busDict
    d[idx] = name
  end
  return d
end

function _branch_result_dict(result::ModelRunResult, branch_labels::Vector{String})::Dict{Tuple{String,String,String},NamedTuple}
  names = _bus_name_by_idx(result.net)
  d = Dict{Tuple{String,String,String},NamedTuple}()
  for row in result.report.branches
    idx = row.branch_index
    from_name = get(names, row.from_bus, string(row.from_bus))
    to_name = get(names, row.to_bus, string(row.to_bus))
    label = idx <= length(branch_labels) ? branch_labels[idx] : String(row.branch)
    nr = _branch_nr_from_label(label)
    d[_branch_key(from_name, to_name, nr)] = merge(row, (label = label, direction = :from, from_name = from_name, to_name = to_name, nr = nr, p_MW = row.p_from_MW, q_MVar = row.q_from_MVar))
    d[_branch_key(to_name, from_name, nr)] = merge(row, (label = label, direction = :to, from_name = to_name, to_name = from_name, nr = nr, p_MW = row.p_to_MW, q_MVar = row.q_to_MVar))
  end
  return d
end

function _compare_bus_rows(ref::For002Scenario, result::ModelRunResult, opt::CliOptions, path::AbstractString)
  calc = _bus_result_dict(result.report)
  rows = NamedTuple[]
  for key in sort(collect(keys(ref.buses)))
    r = ref.buses[key]
    c = get(calc, key, nothing)
    if c === nothing
      push!(rows, (bus_name = r.bus_name, matched = false, d_v_kV = Inf, d_vm_pu = Inf, d_va_deg = Inf, d_p_gen_MW = Inf, d_q_gen_MVar = Inf, d_p_load_MW = Inf, d_q_load_MVar = Inf))
      continue
    end
    ref_vm_pu = r.v_kV / c.vn_kV
    push!(
      rows,
      (
        bus_name = r.bus_name,
        matched = true,
        ref_v_kV = r.v_kV,
        calc_v_kV = c.v_kV,
        d_v_kV = c.v_kV - r.v_kV,
        ref_vm_pu = ref_vm_pu,
        calc_vm_pu = c.vm_pu,
        d_vm_pu = c.vm_pu - ref_vm_pu,
        ref_va_deg = r.va_deg,
        calc_va_deg = c.va_deg,
        d_va_deg = _angle_delta_deg(c.va_deg, r.va_deg),
        ref_p_gen_MW = r.p_gen_MW,
        calc_p_gen_MW = c.p_gen_MW,
        d_p_gen_MW = c.p_gen_MW - r.p_gen_MW,
        ref_q_gen_MVar = r.q_gen_MVar,
        calc_q_gen_MVar = c.q_gen_MVar,
        d_q_gen_MVar = c.q_gen_MVar - r.q_gen_MVar,
        ref_p_load_MW = r.p_load_MW,
        calc_p_load_MW = c.p_load_MW,
        d_p_load_MW = c.p_load_MW - r.p_load_MW,
        ref_q_load_MVar = r.q_load_MVar,
        calc_q_load_MVar = c.q_load_MVar,
        d_q_load_MVar = c.q_load_MVar - r.q_load_MVar,
      ),
    )
  end

  open(path, "w") do io
    println(io, "bus_name,matched,ref_v_kV,calc_v_kV,d_v_kV,ref_vm_pu,calc_vm_pu,d_vm_pu,ref_va_deg,calc_va_deg,d_va_deg,ref_p_gen_MW,calc_p_gen_MW,d_p_gen_MW,ref_q_gen_MVar,calc_q_gen_MVar,d_q_gen_MVar,ref_p_load_MW,calc_p_load_MW,d_p_load_MW,ref_q_load_MVar,calc_q_load_MVar,d_q_load_MVar")
    for row in rows
      if row.matched
        println(
          io,
          _csv_line(
            row.bus_name,
            row.matched,
            row.ref_v_kV,
            row.calc_v_kV,
            row.d_v_kV,
            row.ref_vm_pu,
            row.calc_vm_pu,
            row.d_vm_pu,
            row.ref_va_deg,
            row.calc_va_deg,
            row.d_va_deg,
            row.ref_p_gen_MW,
            row.calc_p_gen_MW,
            row.d_p_gen_MW,
            row.ref_q_gen_MVar,
            row.calc_q_gen_MVar,
            row.d_q_gen_MVar,
            row.ref_p_load_MW,
            row.calc_p_load_MW,
            row.d_p_load_MW,
            row.ref_q_load_MVar,
            row.calc_q_load_MVar,
            row.d_q_load_MVar,
          ),
        )
      else
        println(io, _csv_line(row.bus_name, row.matched, "", "", row.d_v_kV, "", "", row.d_vm_pu, "", "", row.d_va_deg, "", "", row.d_p_gen_MW, "", "", row.d_q_gen_MVar, "", "", row.d_p_load_MW, "", "", row.d_q_load_MVar))
      end
    end
  end

  finite_abs(field) = [abs(getfield(row, field)) for row in rows if isfinite(getfield(row, field))]
  max_or_inf(v) = isempty(v) ? Inf : maximum(v)
  missing_count = count(row -> !row.matched, rows)
  return (
    csv = path,
    count = length(rows),
    missing = missing_count,
    max_abs_d_v_kV = max_or_inf(finite_abs(:d_v_kV)),
    max_abs_d_vm_pu = max_or_inf(finite_abs(:d_vm_pu)),
    max_abs_d_va_deg = max_or_inf(finite_abs(:d_va_deg)),
    max_abs_d_p_gen_MW = max_or_inf(finite_abs(:d_p_gen_MW)),
    max_abs_d_q_gen_MVar = max_or_inf(finite_abs(:d_q_gen_MVar)),
    max_abs_d_p_load_MW = max_or_inf(finite_abs(:d_p_load_MW)),
    max_abs_d_q_load_MVar = max_or_inf(finite_abs(:d_q_load_MVar)),
  )
end

function _compare_branch_rows(ref::For002Scenario, result::ModelRunResult, branch_labels::Vector{String}, opt::CliOptions, path::AbstractString)
  refd = _reference_branch_dict(ref)
  calcd = _branch_result_dict(result, branch_labels)
  keys_all = sort(collect(keys(refd)); by = x -> (x[1], x[2], x[3]))
  rows = NamedTuple[]
  for key in keys_all
    r = refd[key]
    c = get(calcd, key, nothing)
    if c === nothing
      push!(rows, (from_bus = r.from_bus, to_bus = r.to_bus, nr = r.nr, matched = false, d_p_MW = Inf, d_q_MVar = Inf, d_p_loss_MW = Inf, d_q_loss_MVar = Inf))
      continue
    end
    push!(
      rows,
      (
        from_bus = r.from_bus,
        to_bus = r.to_bus,
        nr = r.nr,
        label = c.label,
        matched = true,
        ref_p_MW = r.p_MW,
        calc_p_MW = c.p_MW,
        d_p_MW = c.p_MW - r.p_MW,
        ref_q_MVar = r.q_MVar,
        calc_q_MVar = c.q_MVar,
        d_q_MVar = c.q_MVar - r.q_MVar,
        ref_p_loss_MW = r.p_loss_MW,
        calc_p_loss_MW = c.p_loss_MW,
        d_p_loss_MW = c.p_loss_MW - r.p_loss_MW,
        ref_q_loss_MVar = r.q_loss_MVar,
        calc_q_loss_MVar = c.q_loss_MVar,
        d_q_loss_MVar = c.q_loss_MVar - r.q_loss_MVar,
      ),
    )
  end

  open(path, "w") do io
    println(io, "from_bus,to_bus,nr,label,matched,ref_p_MW,calc_p_MW,d_p_MW,ref_q_MVar,calc_q_MVar,d_q_MVar,ref_p_loss_MW,calc_p_loss_MW,d_p_loss_MW,ref_q_loss_MVar,calc_q_loss_MVar,d_q_loss_MVar")
    for row in rows
      if row.matched
        println(io, _csv_line(row.from_bus, row.to_bus, row.nr, row.label, row.matched, row.ref_p_MW, row.calc_p_MW, row.d_p_MW, row.ref_q_MVar, row.calc_q_MVar, row.d_q_MVar, row.ref_p_loss_MW, row.calc_p_loss_MW, row.d_p_loss_MW, row.ref_q_loss_MVar, row.calc_q_loss_MVar, row.d_q_loss_MVar))
      else
        println(io, _csv_line(row.from_bus, row.to_bus, row.nr, "", row.matched, "", "", row.d_p_MW, "", "", row.d_q_MVar, "", "", row.d_p_loss_MW, "", "", row.d_q_loss_MVar))
      end
    end
  end

  finite_abs(field) = [abs(getfield(row, field)) for row in rows if isfinite(getfield(row, field))]
  max_or_inf(v) = isempty(v) ? Inf : maximum(v)
  return (
    csv = path,
    count = length(rows),
    missing = count(row -> !row.matched, rows),
    max_abs_d_p_MW = max_or_inf(finite_abs(:d_p_MW)),
    max_abs_d_q_MVar = max_or_inf(finite_abs(:d_q_MVar)),
    max_abs_d_p_loss_MW = max_or_inf(finite_abs(:d_p_loss_MW)),
    max_abs_d_q_loss_MVar = max_or_inf(finite_abs(:d_q_loss_MVar)),
  )
end

function _compare_totals(ref::For002Scenario, result::ModelRunResult)
  p_loss = result.report.metadata.total_p_loss_MW
  q_loss = result.report.metadata.total_q_loss_MVar
  d_p = ref.total_p_loss_MW === nothing ? nothing : p_loss - ref.total_p_loss_MW
  d_q = ref.total_q_loss_MVar === nothing ? nothing : q_loss - ref.total_q_loss_MVar
  return (ref_p_loss_MW = ref.total_p_loss_MW, calc_p_loss_MW = p_loss, d_p_loss_MW = d_p, ref_q_loss_MVar = ref.total_q_loss_MVar, calc_q_loss_MVar = q_loss, d_q_loss_MVar = d_q)
end

function _base_pass(bus_cmp, branch_cmp, totals, result::ModelRunResult, opt::CliOptions)::Bool
  total_ok = (totals.d_p_loss_MW === nothing || abs(totals.d_p_loss_MW) <= opt.tol_p_MW) && (totals.d_q_loss_MVar === nothing || abs(totals.d_q_loss_MVar) <= opt.tol_q_MVar)
  return result.status == 0 &&
         bus_cmp.missing == 0 &&
         bus_cmp.max_abs_d_v_kV <= opt.tol_v_kV &&
         bus_cmp.max_abs_d_vm_pu <= opt.tol_vm_pu &&
         bus_cmp.max_abs_d_va_deg <= opt.tol_va_deg &&
         bus_cmp.max_abs_d_p_gen_MW <= opt.tol_p_MW &&
         bus_cmp.max_abs_d_q_gen_MVar <= opt.tol_q_MVar &&
         bus_cmp.max_abs_d_p_load_MW <= opt.tol_p_MW &&
         bus_cmp.max_abs_d_q_load_MVar <= opt.tol_q_MVar &&
         branch_cmp.missing == 0 &&
         branch_cmp.max_abs_d_p_MW <= opt.tol_branch_p_MW &&
         branch_cmp.max_abs_d_q_MVar <= opt.tol_branch_q_MVar &&
         total_ok
end

_diagnostic_objective(bus_cmp, branch_cmp) =
  bus_cmp.max_abs_d_v_kV +
  10.0 * bus_cmp.max_abs_d_va_deg +
  0.05 * bus_cmp.max_abs_d_q_gen_MVar +
  0.05 * branch_cmp.max_abs_d_q_MVar

function _with_branch_b_sweep_point(opt::CliOptions, line_b_scale::Float64, trafo_b_scale::Float64)::CliOptions
  return CliOptions(
    opt.for002_path,
    opt.matpower_path,
    opt.builder_path,
    opt.builder_function,
    opt.output_dir,
    opt.run_matpower,
    opt.run_builder,
    false,
    opt.maxiter,
    opt.solver_tol,
    opt.verbose,
    opt.autodamp,
    opt.opt_flatstart,
    opt.qlimits_enabled,
    opt.tol_v_kV,
    opt.tol_vm_pu,
    opt.tol_va_deg,
    opt.tol_p_MW,
    opt.tol_q_MVar,
    opt.tol_branch_p_MW,
    opt.tol_branch_q_MVar,
    opt.branch_b_scale,
    line_b_scale,
    trafo_b_scale,
    :normal,
    true,
    true,
    opt.sweep_branch_b,
    opt.infer_pv_buses,
    opt.diagnose_active_flow,
    opt.diagnose_worst_branches,
    opt.sweep_worst_branch_angle,
  )
end

function _run_branch_b_sweep(opt::CliOptions, ref_base::For002Scenario, model_builders::Vector{Tuple{String,Function}}, branch_labels::Vector{String})
  line_scales = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
  trafo_scales = [0.0, 0.5, 1.0]
  rows = NamedTuple[]
  scratch_dir = joinpath(opt.output_dir, "_branch_b_sweep_details")
  mkpath(scratch_dir)

  for line_b_scale in line_scales, trafo_b_scale in trafo_scales
    sweep_opt = _with_branch_b_sweep_point(opt, line_b_scale, trafo_b_scale)
    metrics = Dict{String,NamedTuple}()
    for (model_name, buildfn) in model_builders
      result = _run_model(model_name, buildfn, sweep_opt, "base")
      tag = "$(model_name)_line$(_case_tag(string(line_b_scale)))_trafo$(_case_tag(string(trafo_b_scale)))"
      bus_cmp = _compare_bus_rows(ref_base, result, sweep_opt, joinpath(scratch_dir, tag * "_bus_comparison.csv"))
      branch_cmp = _compare_branch_rows(ref_base, result, branch_labels, sweep_opt, joinpath(scratch_dir, tag * "_branch_comparison.csv"))
      totals = _compare_totals(ref_base, result)
      metrics[String(model_name)] = (
        status = result.status,
        max_abs_d_v_kV = bus_cmp.max_abs_d_v_kV,
        max_abs_d_va_deg = bus_cmp.max_abs_d_va_deg,
        max_abs_d_q_gen_MVar = bus_cmp.max_abs_d_q_gen_MVar,
        max_abs_branch_d_p_MW = branch_cmp.max_abs_d_p_MW,
        max_abs_branch_d_q_MVar = branch_cmp.max_abs_d_q_MVar,
        total_d_p_loss_MW = totals.d_p_loss_MW,
        total_d_q_loss_MVar = totals.d_q_loss_MVar,
        objective = _diagnostic_objective(bus_cmp, branch_cmp),
      )
    end
    missing_metrics = (status = missing, max_abs_d_v_kV = missing, max_abs_d_va_deg = missing, max_abs_d_q_gen_MVar = missing, max_abs_branch_d_p_MW = missing, max_abs_branch_d_q_MVar = missing, total_d_p_loss_MW = missing, total_d_q_loss_MVar = missing, objective = Inf)
    matpower = get(metrics, "matpower", missing_metrics)
    builder = get(metrics, "builder", missing_metrics)
    push!(
      rows,
      (
        line_b_scale = line_b_scale,
        trafo_b_scale = trafo_b_scale,
        matpower_status = matpower.status,
        builder_status = builder.status,
        matpower_max_abs_d_v_kV = matpower.max_abs_d_v_kV,
        matpower_max_abs_d_va_deg = matpower.max_abs_d_va_deg,
        matpower_max_abs_d_q_gen_MVar = matpower.max_abs_d_q_gen_MVar,
        matpower_max_abs_branch_d_p_MW = matpower.max_abs_branch_d_p_MW,
        matpower_max_abs_branch_d_q_MVar = matpower.max_abs_branch_d_q_MVar,
        matpower_total_d_p_loss_MW = matpower.total_d_p_loss_MW,
        matpower_total_d_q_loss_MVar = matpower.total_d_q_loss_MVar,
        matpower_objective = matpower.objective,
        builder_max_abs_d_v_kV = builder.max_abs_d_v_kV,
        builder_max_abs_d_va_deg = builder.max_abs_d_va_deg,
        builder_max_abs_d_q_gen_MVar = builder.max_abs_d_q_gen_MVar,
        builder_max_abs_branch_d_p_MW = builder.max_abs_branch_d_p_MW,
        builder_max_abs_branch_d_q_MVar = builder.max_abs_branch_d_q_MVar,
        builder_total_d_p_loss_MW = builder.total_d_p_loss_MW,
        builder_total_d_q_loss_MVar = builder.total_d_q_loss_MVar,
        builder_objective = builder.objective,
      ),
    )
  end
  return rows
end

function _write_branch_b_sweep_summary_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "line_b_scale,trafo_b_scale,matpower_status,builder_status,matpower_max_abs_d_v_kV,matpower_max_abs_d_va_deg,matpower_max_abs_d_q_gen_MVar,matpower_max_abs_branch_d_p_MW,matpower_max_abs_branch_d_q_MVar,matpower_total_d_p_loss_MW,matpower_total_d_q_loss_MVar,matpower_objective,builder_max_abs_d_v_kV,builder_max_abs_d_va_deg,builder_max_abs_d_q_gen_MVar,builder_max_abs_branch_d_p_MW,builder_max_abs_branch_d_q_MVar,builder_total_d_p_loss_MW,builder_total_d_q_loss_MVar,builder_objective")
    for row in rows
      println(io, _csv_line(row.line_b_scale, row.trafo_b_scale, row.matpower_status, row.builder_status, row.matpower_max_abs_d_v_kV, row.matpower_max_abs_d_va_deg, row.matpower_max_abs_d_q_gen_MVar, row.matpower_max_abs_branch_d_p_MW, row.matpower_max_abs_branch_d_q_MVar, row.matpower_total_d_p_loss_MW, row.matpower_total_d_q_loss_MVar, row.matpower_objective, row.builder_max_abs_d_v_kV, row.builder_max_abs_d_va_deg, row.builder_max_abs_d_q_gen_MVar, row.builder_max_abs_branch_d_p_MW, row.builder_max_abs_branch_d_q_MVar, row.builder_total_d_p_loss_MW, row.builder_total_d_q_loss_MVar, row.builder_objective))
    end
  end
  return path
end

function _print_best_sweep_rows(rows::Vector{NamedTuple})
  isempty(rows) && return nothing
  println("Branch-b sweep objective: max_abs_d_v_kV + 10.0*max_abs_d_va_deg + 0.05*max_abs_d_q_gen_MVar + 0.05*max_abs_branch_d_q_MVar")
  for (label, field) in (("MATPOWER", :matpower_objective), ("Builder", :builder_objective))
    best = first(sort(rows; by = row -> getfield(row, field)))
    @printf(
      "  Best %-8s line_b_scale=%.2f trafo_b_scale=%.1f objective=%.6g max|dV|=%.6g kV max|dVa|=%.6g deg max|dQgen|=%.6g MVar max branch |dQ|=%.6g MVar\n",
      label,
      best.line_b_scale,
      best.trafo_b_scale,
      getfield(best, field),
      getfield(best, Symbol(lowercase(label) * "_max_abs_d_v_kV")),
      getfield(best, Symbol(lowercase(label) * "_max_abs_d_va_deg")),
      getfield(best, Symbol(lowercase(label) * "_max_abs_d_q_gen_MVar")),
      getfield(best, Symbol(lowercase(label) * "_max_abs_branch_d_q_MVar")),
    )
  end
  println()
  return nothing
end

function _compare_two_model_runs(a::ModelRunResult, b::ModelRunResult, path::AbstractString)
  ad = _bus_result_dict(a.report)
  bd = _bus_result_dict(b.report)
  keys_all = sort(collect(union(keys(ad), keys(bd))))
  open(path, "w") do io
    println(io, "bus_name,model_a,model_b,matched,vm_a_pu,vm_b_pu,d_vm_pu,va_a_deg,va_b_deg,d_va_deg,v_a_kV,v_b_kV,d_v_kV,p_gen_a_MW,p_gen_b_MW,d_p_gen_MW,q_gen_a_MVar,q_gen_b_MVar,d_q_gen_MVar")
    for key in keys_all
      ra = get(ad, key, nothing)
      rb = get(bd, key, nothing)
      if ra === nothing || rb === nothing
        println(io, _csv_line(key, a.model_name, b.model_name, false, "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""))
        continue
      end
      println(
        io,
        _csv_line(
          ra.bus_name,
          a.model_name,
          b.model_name,
          true,
          ra.vm_pu,
          rb.vm_pu,
          rb.vm_pu - ra.vm_pu,
          ra.va_deg,
          rb.va_deg,
          _angle_delta_deg(rb.va_deg, ra.va_deg),
          ra.v_kV,
          rb.v_kV,
          rb.v_kV - ra.v_kV,
          ra.p_gen_MW,
          rb.p_gen_MW,
          rb.p_gen_MW - ra.p_gen_MW,
          ra.q_gen_MVar,
          rb.q_gen_MVar,
          rb.q_gen_MVar - ra.q_gen_MVar,
        ),
      )
    end
  end
  return path
end

function _worst_offender_rows(ref::For002Scenario, result::ModelRunResult, branch_labels::Vector{String}; top_n::Int = 40)
  rows = NamedTuple[]
  calc_buses = _bus_result_dict(result.report)
  for key in sort(collect(keys(ref.buses)))
    r = ref.buses[key]
    c = get(calc_buses, key, nothing)
    c === nothing && continue
    bus_name = String(r.bus_name)
    candidates = (
      ("d_v_kV", c.v_kV - r.v_kV, r.v_kV, c.v_kV),
      ("d_va_deg", _angle_delta_deg(c.va_deg, r.va_deg), r.va_deg, c.va_deg),
      ("d_p_gen_MW", c.p_gen_MW - r.p_gen_MW, r.p_gen_MW, c.p_gen_MW),
      ("d_q_gen_MVar", c.q_gen_MVar - r.q_gen_MVar, r.q_gen_MVar, c.q_gen_MVar),
      ("d_p_load_MW", c.p_load_MW - r.p_load_MW, r.p_load_MW, c.p_load_MW),
      ("d_q_load_MVar", c.q_load_MVar - r.q_load_MVar, r.q_load_MVar, c.q_load_MVar),
    )
    for (deviation_type, signed_deviation, reference_value, calculated_value) in candidates
      push!(rows, (kind = "bus", model = result.model_name, element = bus_name, deviation_type = deviation_type, absolute_deviation = abs(signed_deviation), signed_deviation = signed_deviation, reference_value = reference_value, calculated_value = calculated_value))
    end
  end

  ref_branches = _reference_branch_dict(ref)
  calc_branches = _branch_result_dict(result, branch_labels)
  for key in sort(collect(keys(ref_branches)); by = x -> (x[1], x[2], x[3]))
    r = ref_branches[key]
    c = get(calc_branches, key, nothing)
    c === nothing && continue
    element = isempty(c.nr) ? "$(c.from_name) -> $(c.to_name)" : "$(c.from_name) -> $(c.to_name) [$(c.nr)]"
    candidates = (
      ("d_p_MW", c.p_MW - r.p_MW, r.p_MW, c.p_MW),
      ("d_q_MVar", c.q_MVar - r.q_MVar, r.q_MVar, c.q_MVar),
      ("d_p_loss_MW", c.p_loss_MW - r.p_loss_MW, r.p_loss_MW, c.p_loss_MW),
      ("d_q_loss_MVar", c.q_loss_MVar - r.q_loss_MVar, r.q_loss_MVar, c.q_loss_MVar),
    )
    for (deviation_type, signed_deviation, reference_value, calculated_value) in candidates
      push!(rows, (kind = "branch", model = result.model_name, element = element, deviation_type = deviation_type, absolute_deviation = abs(signed_deviation), signed_deviation = signed_deviation, reference_value = reference_value, calculated_value = calculated_value))
    end
  end
  sorted_rows = sort(rows; by = row -> row.absolute_deviation, rev = true)
  return sorted_rows[1:min(top_n, length(sorted_rows))]
end

function _write_worst_offenders_csv(path::AbstractString, ref::For002Scenario, base_results::Vector{ModelRunResult}, branch_labels::Vector{String})
  rows = NamedTuple[]
  for result in base_results
    append!(rows, _worst_offender_rows(ref, result, branch_labels))
  end
  rows = sort(rows; by = row -> row.absolute_deviation, rev = true)
  open(path, "w") do io
    println(io, "kind,model,element,deviation_type,absolute_deviation,signed_deviation,reference_value,calculated_value")
    for row in rows
      println(io, _csv_line(row.kind, row.model, row.element, row.deviation_type, row.absolute_deviation, row.signed_deviation, row.reference_value, row.calculated_value))
    end
  end
  return path
end

function _matpower_generator_bus_metadata(mpc)
  by_bus = Dict{Int,NamedTuple}()
  size(mpc.gen, 2) >= 10 || return by_bus
  @inbounds for g in axes(mpc.gen, 1)
    size(mpc.gen, 2) >= 8 && Float64(mpc.gen[g, 8]) <= 0.0 && continue
    bus = Int(mpc.gen[g, 1])
    prev = get(by_bus, bus, (qmin = nothing, qmax = nothing, vg = nothing))
    qmin = Float64(mpc.gen[g, 5])
    qmax = Float64(mpc.gen[g, 4])
    vg = size(mpc.gen, 2) >= 6 ? Float64(mpc.gen[g, 6]) : NaN
    by_bus[bus] = (
      qmin = prev.qmin === nothing ? qmin : prev.qmin + qmin,
      qmax = prev.qmax === nothing ? qmax : prev.qmax + qmax,
      vg = prev.vg === nothing && isfinite(vg) && vg > 0.0 ? vg : prev.vg,
    )
  end
  return by_bus
end

function _first_valid_generator_vset(generators)::Union{Nothing,Float64}
  for ps in generators
    ps.vm_pu === nothing && continue
    vm = Float64(ps.vm_pu)
    isfinite(vm) && vm > 0.0 && return vm
  end
  return nothing
end

function _pv_inference_rows_and_apply!(net::Net, opt::CliOptions, model_name::AbstractString; mpc = nothing, apply::Bool = true)
  mode = opt.infer_pv_buses
  rows = NamedTuple[]
  bus_names = _bus_name_by_idx(net)
  matpower_meta = mpc === nothing ? Dict{Int,NamedTuple}() : _matpower_generator_bus_metadata(mpc)

  for bus_idx in sort(collect(keys(bus_names)))
    node = net.nodeVec[bus_idx]
    prosumers = Sparlectra.getBusProsumers(net, bus_idx)
    generators = [ps for ps in prosumers if Sparlectra.isGenerator(ps)]
    isempty(generators) && continue

    original_type = String(Sparlectra.toString(Sparlectra.getNodeType(node)))
    original_regulated = any(ps -> ps.isRegulated, generators)
    p_gen = sum((ps.pVal === nothing ? 0.0 : Float64(ps.pVal)) for ps in generators)
    q_gen = sum((ps.qVal === nothing ? 0.0 : Float64(ps.qVal)) for ps in generators)
    q_min_local = all(ps -> ps.minQ !== nothing, generators) ? sum(Float64(ps.minQ) for ps in generators) : nothing
    q_max_local = all(ps -> ps.maxQ !== nothing, generators) ? sum(Float64(ps.maxQ) for ps in generators) : nothing
    orig_bus = get(net.busOrigIdxDict, bus_idx, bus_idx)
    raw_meta = get(matpower_meta, orig_bus, nothing)
    q_min = raw_meta === nothing ? q_min_local : raw_meta.qmin
    q_max = raw_meta === nothing ? q_max_local : raw_meta.qmax
    vg = raw_meta === nothing ? nothing : raw_meta.vg
    ps_vset = _first_valid_generator_vset(generators)
    current_vm = isfinite(node._vm_pu) && node._vm_pu > 0.0 ? node._vm_pu : nothing

    vset = nothing
    should_infer = false
    reason = "mode off"
    if mode == :off
      vset = vg === nothing ? ps_vset : vg
    elseif Sparlectra.getNodeType(node) == Sparlectra.Slack
      vset = vg === nothing ? something(ps_vset, current_vm) : vg
      reason = "slack generator bus left unchanged"
    elseif mode == Symbol("from-q-limits")
      vset = vg === nothing ? ps_vset : vg
      finite_limits = q_min !== nothing && q_max !== nothing && isfinite(q_min) && isfinite(q_max)
      should_infer = finite_limits && vset !== nothing && isfinite(vset) && vset > 0.0
      reason = should_infer ? "finite Q limits and voltage setpoint" : "missing finite Q limits or voltage setpoint"
    elseif mode == Symbol("from-gen-vg")
      vset = vg === nothing ? ps_vset : vg
      should_infer = vset !== nothing && isfinite(vset) && vset > 0.0
      reason = raw_meta === nothing ? "native generator voltage setpoint" : "MATPOWER GEN.VG voltage setpoint"
    elseif mode == Symbol("all-generator-buses")
      vset = something(vg, ps_vset, current_vm, 1.0)
      should_infer = isfinite(vset) && vset > 0.0
      reason = "diagnostic all non-slack generator buses"
    end

    if apply && should_infer
      for ps in generators
        ps.isRegulated = true
        ps.vm_pu = vset
      end
      Sparlectra.setVmVa!(node = node, vm_pu = vset)
    end

    inferred_type = if should_infer
      "PV"
    else
      original_type
    end
    push!(
      rows,
      (
        model = String(model_name),
        bus_name = bus_names[bus_idx],
        original_bus_type = original_type,
        original_regulation_state = original_regulated,
        inferred_bus_type = inferred_type,
        has_generator = true,
        p_gen_MW = p_gen,
        q_gen_MVar_before = q_gen,
        q_min_MVar = q_min,
        q_max_MVar = q_max,
        voltage_setpoint_pu = vset,
        inference_reason = reason,
      ),
    )
  end

  if apply && mode != :off
    Sparlectra.refreshBusTypesFromProsumers!(net)
    Sparlectra.buildQLimits!(net)
  end
  return rows
end

function _write_pv_bus_diagnostics_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,bus_name,original_bus_type,original_regulation_state,inferred_bus_type,has_generator,p_gen_MW,q_gen_MVar_before,q_min_MVar,q_max_MVar,voltage_setpoint_pu,inference_reason")
    for row in rows
      println(io, _csv_line(row.model, row.bus_name, row.original_bus_type, row.original_regulation_state, row.inferred_bus_type, row.has_generator, row.p_gen_MW, row.q_gen_MVar_before, row.q_min_MVar, row.q_max_MVar, row.voltage_setpoint_pu, row.inference_reason))
    end
  end
  return path
end

function _pv_bus_diagnostics_note(opt::CliOptions, rows::Vector{NamedTuple})::Union{Nothing,String}
  opt.infer_pv_buses == :off && return nothing
  builder_rows = [row for row in rows if row.model == "builder"]
  isempty(builder_rows) && return "Native Builder PV inference was not evaluated because the Builder model path was disabled."
  blocked = [row.bus_name for row in builder_rows if occursin("missing finite Q limits", row.inference_reason)]
  if !isempty(blocked)
    return "Native Builder from-q-limits inference was not technically feasible for buses without finite builder-side Q limits: " * join(blocked, ", ")
  end
  return nothing
end

function _bus_angle_deviation_rows(ref::For002Scenario, result::ModelRunResult)
  calc = _bus_result_dict(result.report)
  rows = NamedTuple[]
  for key in sort(collect(keys(ref.buses)))
    r = ref.buses[key]
    c = get(calc, key, nothing)
    c === nothing && continue
    push!(rows, (bus_name = r.bus_name, d_va_deg = _angle_delta_deg(c.va_deg, r.va_deg), ref_va_deg = r.va_deg, calc_va_deg = c.va_deg))
  end
  return rows
end

function _branch_flow_deviation_rows(ref::For002Scenario, result::ModelRunResult, branch_labels::Vector{String})
  refd = _reference_branch_dict(ref)
  calcd = _branch_result_dict(result, branch_labels)
  rows = NamedTuple[]
  for key in sort(collect(keys(refd)); by = x -> (x[1], x[2], x[3]))
    r = refd[key]
    c = get(calcd, key, nothing)
    c === nothing && continue
    push!(
      rows,
      (
        branch_index = c.branch_index,
        label = c.label,
        from_bus = c.from_name,
        to_bus = c.to_name,
        nr = c.nr,
        d_p_MW = c.p_MW - r.p_MW,
        d_q_MVar = c.q_MVar - r.q_MVar,
        ref_p_MW = r.p_MW,
        calc_p_MW = c.p_MW,
        ref_q_MVar = r.q_MVar,
        calc_q_MVar = c.q_MVar,
      ),
    )
  end
  return rows
end

function _slack_bus_power_deviation(ref::For002Scenario, result::ModelRunResult)
  calc = _bus_result_dict(result.report)
  for node in result.net.nodeVec
    Sparlectra.getNodeType(node) == Sparlectra.Slack || continue
    bus_name = get(_bus_name_by_idx(result.net), node.busIdx, string(node.busIdx))
    r = get(ref.buses, _norm_name(bus_name), nothing)
    c = get(calc, _norm_name(bus_name), nothing)
    r === nothing || c === nothing || return (
      bus_name = bus_name,
      d_p_gen_MW = c.p_gen_MW - r.p_gen_MW,
      d_q_gen_MVar = c.q_gen_MVar - r.q_gen_MVar,
      ref_p_gen_MW = r.p_gen_MW,
      calc_p_gen_MW = c.p_gen_MW,
      ref_q_gen_MVar = r.q_gen_MVar,
      calc_q_gen_MVar = c.q_gen_MVar,
    )
  end
  return nothing
end

function _active_power_balance_rows(ref::For002Scenario, result::ModelRunResult)
  rows = NamedTuple[]
  ref_buses = ref.buses
  for row in result.report.nodes
    key = _norm_name(row.bus_name)
    refrow = get(ref_buses, key, nothing)
    node_idx = get(result.net.busDict, row.bus_name, 0)
    bus_type = node_idx == 0 ? "" : String(Sparlectra.toString(Sparlectra.getNodeType(result.net.nodeVec[node_idx])))
    bus_prosumers = node_idx == 0 ? Any[] : Sparlectra.getBusProsumers(result.net, node_idx)
    regulated = any(ps -> ps.isRegulated, bus_prosumers)
    is_slack_bus = node_idx != 0 && Sparlectra.getNodeType(result.net.nodeVec[node_idx]) == Sparlectra.Slack
    push!(
      rows,
      (
        model = result.model_name,
        bus_name = row.bus_name,
        bus_type = bus_type,
        regulation_state = regulated,
        p_gen_MW = row.p_gen_MW,
        q_gen_MVar = row.q_gen_MVar,
        p_load_MW = row.p_load_MW,
        q_load_MVar = row.q_load_MVar,
        net_p_injection_MW = row.p_gen_MW - row.p_load_MW,
        net_q_injection_MVar = row.q_gen_MVar - row.q_load_MVar,
        ref_p_gen_MW = refrow === nothing ? nothing : refrow.p_gen_MW,
        ref_q_gen_MVar = refrow === nothing ? nothing : refrow.q_gen_MVar,
        ref_p_load_MW = refrow === nothing ? nothing : refrow.p_load_MW,
        ref_q_load_MVar = refrow === nothing ? nothing : refrow.q_load_MVar,
        d_p_gen_MW = refrow === nothing ? nothing : row.p_gen_MW - refrow.p_gen_MW,
        d_q_gen_MVar = refrow === nothing ? nothing : row.q_gen_MVar - refrow.q_gen_MVar,
        d_p_load_MW = refrow === nothing ? nothing : row.p_load_MW - refrow.p_load_MW,
        d_q_load_MVar = refrow === nothing ? nothing : row.q_load_MVar - refrow.q_load_MVar,
        is_slack_bus = is_slack_bus,
      ),
    )
  end
  return rows
end

function _voltage_coupling_note(from_kv, to_kv)::String
  if from_kv === missing || to_kv === missing
    return "unknown voltage coupling"
  end
  if isapprox(max(from_kv, to_kv), 400.0; atol = 1.0) && isapprox(min(from_kv, to_kv), 230.0; atol = 1.0)
    return "400/230-kV coupling"
  elseif isapprox(from_kv, to_kv; atol = 1e-6)
    return "same-voltage transformer"
  else
    return "other voltage coupling"
  end
end

function _branch_parameter_rows(result::ModelRunResult, branch_labels::Vector{String}, top_p_indices::Set{Int}, top_angle_branch_indices::Set{Int})
  names = _bus_name_by_idx(result.net)
  rows = NamedTuple[]
  for br in result.net.branchVec
    idx = br.branchIdx
    from_base_kv = result.net.nodeVec[br.fromBus].comp.cVN
    to_base_kv = result.net.nodeVec[br.toBus].comp.cVN
    push!(
      rows,
      (
        model = result.model_name,
        branch_index = idx,
        branch_label = idx <= length(branch_labels) ? branch_labels[idx] : string(idx),
        from_bus = get(names, br.fromBus, string(br.fromBus)),
        to_bus = get(names, br.toBus, string(br.toBus)),
        branch_kind = _branch_model_kind_label(result.net, br),
        r_pu = br.r_pu,
        x_pu = br.x_pu,
        b_pu_after_effective_scaling = br.b_pu,
        ratio = br.ratio,
        angle = br.angle,
        status = br.status,
        from_base_kV = from_base_kv,
        to_base_kV = to_base_kv,
        is_top_p_flow_offender = idx in top_p_indices,
        is_top_angle_offender = idx in top_angle_branch_indices,
      ),
    )
  end
  return rows
end

function _transformer_tap_rows(result::ModelRunResult, branch_labels::Vector{String}, branch_deviations::Vector{NamedTuple})
  names = _bus_name_by_idx(result.net)
  deviation_by_index = Dict(row.branch_index => row for row in branch_deviations)
  rows = NamedTuple[]
  for br in result.net.branchVec
    _branch_diagnostic_kind(result.net, br) == :transformer || continue
    idx = br.branchIdx
    from_base_kv = result.net.nodeVec[br.fromBus].comp.cVN
    to_base_kv = result.net.nodeVec[br.toBus].comp.cVN
    dev = get(deviation_by_index, idx, nothing)
    push!(
      rows,
      (
        model = result.model_name,
        branch_index = idx,
        branch_label = idx <= length(branch_labels) ? branch_labels[idx] : string(idx),
        from_bus = get(names, br.fromBus, string(br.fromBus)),
        to_bus = get(names, br.toBus, string(br.toBus)),
        from_base_kV = from_base_kv,
        to_base_kV = to_base_kv,
        ratio = br.ratio,
        angle_deg = br.angle,
        r_pu = br.r_pu,
        x_pu = br.x_pu,
        b_pu = br.b_pu,
        d_p_MW = dev === nothing ? nothing : dev.d_p_MW,
        d_q_MVar = dev === nothing ? nothing : dev.d_q_MVar,
        coupling_note = _voltage_coupling_note(from_base_kv, to_base_kv),
      ),
    )
  end
  return rows
end

function _write_active_power_balance_diagnostics_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,bus_name,bus_type,regulation_state,p_gen_MW,q_gen_MVar,p_load_MW,q_load_MVar,net_p_injection_MW,net_q_injection_MVar,ref_p_gen_MW,ref_q_gen_MVar,ref_p_load_MW,ref_q_load_MVar,d_p_gen_MW,d_q_gen_MVar,d_p_load_MW,d_q_load_MVar,is_slack_bus")
    for row in rows
      println(io, _csv_line(row.model, row.bus_name, row.bus_type, row.regulation_state, row.p_gen_MW, row.q_gen_MVar, row.p_load_MW, row.q_load_MVar, row.net_p_injection_MW, row.net_q_injection_MVar, row.ref_p_gen_MW, row.ref_q_gen_MVar, row.ref_p_load_MW, row.ref_q_load_MVar, row.d_p_gen_MW, row.d_q_gen_MVar, row.d_p_load_MW, row.d_q_load_MVar, row.is_slack_bus))
    end
  end
  return path
end

function _write_branch_parameter_diagnostics_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,branch_index,branch_label,from_bus,to_bus,branch_kind,r_pu,x_pu,b_pu_after_effective_scaling,ratio,angle,status,from_base_kV,to_base_kV,is_top_p_flow_offender,is_top_angle_offender")
    for row in rows
      println(io, _csv_line(row.model, row.branch_index, row.branch_label, row.from_bus, row.to_bus, row.branch_kind, row.r_pu, row.x_pu, row.b_pu_after_effective_scaling, row.ratio, row.angle, row.status, row.from_base_kV, row.to_base_kV, row.is_top_p_flow_offender, row.is_top_angle_offender))
    end
  end
  return path
end

function _write_transformer_tap_diagnostics_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,branch_index,branch_label,from_bus,to_bus,from_base_kV,to_base_kV,ratio,angle_deg,r_pu,x_pu,b_pu,d_p_MW,d_q_MVar,coupling_note")
    for row in rows
      println(io, _csv_line(row.model, row.branch_index, row.branch_label, row.from_bus, row.to_bus, row.from_base_kV, row.to_base_kV, row.ratio, row.angle_deg, row.r_pu, row.x_pu, row.b_pu, row.d_p_MW, row.d_q_MVar, row.coupling_note))
    end
  end
  return path
end

function _write_active_flow_diagnostics(opt::CliOptions, ref::For002Scenario, base_results::Vector{ModelRunResult}, branch_labels::Vector{String})
  active_rows = NamedTuple[]
  branch_parameter_rows = NamedTuple[]
  transformer_rows = NamedTuple[]
  console_rows = NamedTuple[]
  for result in base_results
    branch_deviations = _branch_flow_deviation_rows(ref, result, branch_labels)
    angle_deviations = _bus_angle_deviation_rows(ref, result)
    top_p = first(sort(branch_deviations; by = row -> abs(row.d_p_MW), rev = true), min(10, length(branch_deviations)))
    top_angles = first(sort(angle_deviations; by = row -> abs(row.d_va_deg), rev = true), min(10, length(angle_deviations)))
    top_p_indices = Set(row.branch_index for row in top_p)
    top_angle_buses = Set(_norm_name(row.bus_name) for row in top_angles)
    top_angle_branch_indices = Set(
      br.branchIdx for br in result.net.branchVec if
      _norm_name(get(_bus_name_by_idx(result.net), br.fromBus, string(br.fromBus))) in top_angle_buses ||
      _norm_name(get(_bus_name_by_idx(result.net), br.toBus, string(br.toBus))) in top_angle_buses
    )
    append!(active_rows, _active_power_balance_rows(ref, result))
    append!(branch_parameter_rows, _branch_parameter_rows(result, branch_labels, top_p_indices, top_angle_branch_indices))
    append!(transformer_rows, _transformer_tap_rows(result, branch_labels, branch_deviations))
    push!(console_rows, (model = result.model_name, top_p = top_p, top_angles = top_angles, slack = _slack_bus_power_deviation(ref, result)))
  end
  active_csv = _write_active_power_balance_diagnostics_csv(joinpath(opt.output_dir, "active_power_balance_diagnostics.csv"), active_rows)
  branch_csv = _write_branch_parameter_diagnostics_csv(joinpath(opt.output_dir, "branch_parameter_diagnostics.csv"), branch_parameter_rows)
  transformer_csv = _write_transformer_tap_diagnostics_csv(joinpath(opt.output_dir, "transformer_tap_diagnostics.csv"), transformer_rows)
  return (active_csv = active_csv, branch_csv = branch_csv, transformer_csv = transformer_csv, console_rows = console_rows)
end

function _branch_by_index(net::Net, idx::Int)
  for br in net.branchVec
    br.branchIdx == idx && return br
  end
  return nothing
end

function _estimated_dc_p_mw(net::Net, br, from_idx::Int, to_idx::Int)
  (br === nothing || !isfinite(br.x_pu) || abs(br.x_pu) < eps(Float64)) && return nothing
  va_from = net.nodeVec[from_idx]._va_deg
  va_to = net.nodeVec[to_idx]._va_deg
  theta = deg2rad(_angle_delta_deg(va_from - br.angle, va_to))
  return net.baseMVA * theta / br.x_pu
end

function _worst_branch_diagnostic_rows(ref::For002Scenario, result::ModelRunResult, branch_labels::Vector{String})
  deviations = _branch_flow_deviation_rows(ref, result, branch_labels)
  top = first(sort(deviations; by = row -> abs(row.d_p_MW), rev = true), min(10, length(deviations)))
  top_keys = Set((row.branch_index, _norm_name(row.from_bus), _norm_name(row.to_bus)) for row in top)
  refd = _reference_branch_dict(ref)
  calcd = _branch_result_dict(result, branch_labels)
  rows = NamedTuple[]
  for key in sort(collect(keys(refd)); by = x -> (x[1], x[2], x[3]))
    r = refd[key]
    c = get(calcd, key, nothing)
    c === nothing && continue
    br = _branch_by_index(result.net, c.branch_index)
    from_idx = c.direction == :from ? c.from_bus : c.to_bus
    to_idx = c.direction == :from ? c.to_bus : c.from_bus
    from_node = result.net.nodeVec[from_idx]
    to_node = result.net.nodeVec[to_idx]
    from_va = from_node._va_deg
    to_va = to_node._va_deg
    push!(
      rows,
      (
        model = result.model_name,
        branch_index = c.branch_index,
        branch_label = c.label,
        comparison_direction_from = c.from_name,
        comparison_direction_to = c.to_name,
        branch_kind = br === nothing ? "missing" : _branch_model_kind_label(result.net, br),
        from_bus_internal = from_idx,
        to_bus_internal = to_idx,
        from_base_kV = from_node.comp.cVN,
        to_base_kV = to_node.comp.cVN,
        r_pu = br === nothing ? missing : br.r_pu,
        x_pu = br === nothing ? missing : br.x_pu,
        b_pu_after_scaling = br === nothing ? missing : br.b_pu,
        ratio = br === nothing ? missing : br.ratio,
        angle_deg = br === nothing ? missing : br.angle,
        status = br === nothing ? missing : br.status,
        calculated_p_MW = c.p_MW,
        reference_for002_p_MW = r.p_MW,
        d_p_MW = c.p_MW - r.p_MW,
        calculated_q_MVar = c.q_MVar,
        reference_for002_q_MVar = r.q_MVar,
        d_q_MVar = c.q_MVar - r.q_MVar,
        calculated_loss_p_MW = c.p_loss_MW,
        reference_loss_p_MW = r.p_loss_MW,
        calculated_loss_q_MVar = c.q_loss_MVar,
        reference_loss_q_MVar = r.q_loss_MVar,
        from_va_deg = from_va,
        to_va_deg = to_va,
        angle_difference_deg = _angle_delta_deg(from_va, to_va),
        estimated_dc_p_MW = br === nothing ? nothing : _estimated_dc_p_mw(result.net, br, from_idx, to_idx),
        is_top_active_flow_offender = (c.branch_index, _norm_name(c.from_name), _norm_name(c.to_name)) in top_keys,
      ),
    )
  end
  return sort(rows; by = row -> abs(row.d_p_MW), rev = true)
end

function _write_worst_branch_flow_diagnostics_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,branch_index,branch_label,comparison_direction_from,comparison_direction_to,branch_kind,from_bus_internal,to_bus_internal,from_base_kV,to_base_kV,r_pu,x_pu,b_pu_after_scaling,ratio,angle_deg,status,calculated_p_MW,reference_for002_p_MW,d_p_MW,calculated_q_MVar,reference_for002_q_MVar,d_q_MVar,calculated_loss_p_MW,reference_loss_p_MW,calculated_loss_q_MVar,reference_loss_q_MVar,from_va_deg,to_va_deg,angle_difference_deg,estimated_dc_p_MW,is_top_active_flow_offender")
    for row in rows
      println(io, _csv_line(row.model, row.branch_index, row.branch_label, row.comparison_direction_from, row.comparison_direction_to, row.branch_kind, row.from_bus_internal, row.to_bus_internal, row.from_base_kV, row.to_base_kV, row.r_pu, row.x_pu, row.b_pu_after_scaling, row.ratio, row.angle_deg, row.status, row.calculated_p_MW, row.reference_for002_p_MW, row.d_p_MW, row.calculated_q_MVar, row.reference_for002_q_MVar, row.d_q_MVar, row.calculated_loss_p_MW, row.reference_loss_p_MW, row.calculated_loss_q_MVar, row.reference_loss_q_MVar, row.from_va_deg, row.to_va_deg, row.angle_difference_deg, row.estimated_dc_p_MW, row.is_top_active_flow_offender))
    end
  end
  return path
end

function _run_worst_branch_angle_sweep(opt::CliOptions, ref::For002Scenario, model_builders::Vector{Tuple{String,Function}}, branch_labels::Vector{String}, diagnostic_rows::Vector{NamedTuple})
  offsets = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
  top_indices = unique(row.branch_index for row in diagnostic_rows if row.is_top_active_flow_offender)
  rows = NamedTuple[]
  scratch = joinpath(opt.output_dir, "_worst_branch_angle_sweep_details")
  mkpath(scratch)
  for (model_name, buildfn) in model_builders, idx in top_indices, offset in offsets
    net = buildfn()
    _scale_branch_b!(net, opt)
    br = _branch_by_index(net, idx)
    br !== nothing && (br.angle += offset)
    report, iters, status, residual = _run_powerflow!(net, opt)
    result = ModelRunResult(String(model_name), "base", nothing, net, report, iters, status, residual)
    tag = "$(model_name)_branch$(idx)_angle$(_case_tag(string(offset)))"
    bus_cmp = _compare_bus_rows(ref, result, opt, joinpath(scratch, tag * "_bus_comparison.csv"))
    branch_cmp = _compare_branch_rows(ref, result, branch_labels, opt, joinpath(scratch, tag * "_branch_comparison.csv"))
    totals = _compare_totals(ref, result)
    push!(rows, (model = String(model_name), branch_index = idx, branch_label = idx <= length(branch_labels) ? branch_labels[idx] : string(idx), angle_offset_deg = offset, solver_status = status, max_abs_branch_d_p_MW = branch_cmp.max_abs_d_p_MW, max_abs_d_va_deg = bus_cmp.max_abs_d_va_deg, max_abs_d_v_kV = bus_cmp.max_abs_d_v_kV, total_d_p_loss_MW = totals.d_p_loss_MW, objective = bus_cmp.max_abs_d_v_kV + 10.0 * bus_cmp.max_abs_d_va_deg + branch_cmp.max_abs_d_p_MW))
  end
  return rows
end

function _write_worst_branch_angle_sweep_csv(path::AbstractString, rows::Vector{NamedTuple})
  open(path, "w") do io
    println(io, "model,branch_index,branch_label,angle_offset_deg,solver_status,max_abs_branch_d_p_MW,max_abs_d_va_deg,max_abs_d_v_kV,total_d_p_loss_MW,objective")
    for row in rows
      println(io, _csv_line(row.model, row.branch_index, row.branch_label, row.angle_offset_deg, row.solver_status, row.max_abs_branch_d_p_MW, row.max_abs_d_va_deg, row.max_abs_d_v_kV, row.total_d_p_loss_MW, row.objective))
    end
  end
  return path
end

function _write_worst_branch_diagnostics(opt::CliOptions, ref::For002Scenario, base_results::Vector{ModelRunResult}, model_builders::Vector{Tuple{String,Function}}, branch_labels::Vector{String})
  rows = NamedTuple[]
  for result in base_results
    append!(rows, _worst_branch_diagnostic_rows(ref, result, branch_labels))
  end
  rows = sort(rows; by = row -> abs(row.d_p_MW), rev = true)
  flow_csv = _write_worst_branch_flow_diagnostics_csv(joinpath(opt.output_dir, "worst_branch_flow_diagnostics.csv"), rows)
  sweep_rows = opt.sweep_worst_branch_angle ? _run_worst_branch_angle_sweep(opt, ref, model_builders, branch_labels, rows) : NamedTuple[]
  sweep_csv = opt.sweep_worst_branch_angle ? _write_worst_branch_angle_sweep_csv(joinpath(opt.output_dir, "worst_branch_angle_sweep.csv"), sweep_rows) : nothing
  return (flow_csv = flow_csv, sweep_csv = sweep_csv, rows = rows, sweep_rows = sweep_rows)
end

function _print_worst_branch_console_summary(diagnostics)
  diagnostics === nothing && return nothing
  println("Worst branch-flow diagnostic summary")
  top = first(diagnostics.rows, min(10, length(diagnostics.rows)))
  kind_counts = Dict{String,Int}()
  for row in top
    kind_counts[row.branch_kind] = get(kind_counts, row.branch_kind, 0) + 1
    @printf("  %-8s #%-3d %-34s %-12s -> %-12s dP=% .6g MW r=% .6g x=% .6g ratio=%s angle=%s kind=%s\n", row.model, row.branch_index, row.branch_label, row.comparison_direction_from, row.comparison_direction_to, row.d_p_MW, row.r_pu, row.x_pu, string(row.ratio), string(row.angle_deg), row.branch_kind)
  end
  println("  Top-offender branch kinds: ", join(["$k=$v" for (k, v) in sort(collect(kind_counts))], ", "))
  tapped = [row for row in top if row.ratio !== missing && (abs(Float64(row.ratio)) > 0.0 || abs(Float64(row.angle_deg)) > 0.0)]
  println("  Tapped/phase-shifting among top offenders: ", isempty(tapped) ? "none" : join(unique(row.branch_label for row in tapped), "; "))
  if diagnostics.sweep_csv !== nothing && !isempty(diagnostics.sweep_rows)
    best = first(sort(diagnostics.sweep_rows; by = row -> row.objective))
    @printf("  Best angle sweep: model=%s branch=%d %s offset=% .3g deg max|dP|=% .6g MW max|dVa|=% .6g deg objective=% .6g\n", best.model, best.branch_index, best.branch_label, best.angle_offset_deg, best.max_abs_branch_d_p_MW, best.max_abs_d_va_deg, best.objective)
  end
  println()
  return nothing
end

function _print_active_flow_console_summary(diagnostics)
  diagnostics === nothing && return nothing
  println("Active-flow diagnostic summary")
  for row in diagnostics.console_rows
    println("  Model: ", row.model)
    println("    Top branch |dP|:")
    for br in row.top_p
      @printf("      %-34s %s -> %s dP=% .6g MW dQ=% .6g MVar\n", br.label, br.from_bus, br.to_bus, br.d_p_MW, br.d_q_MVar)
    end
    println("    Top bus |dVa|:")
    for bus in row.top_angles
      @printf("      %-12s dVa=% .6g deg ref=% .6g calc=% .6g\n", bus.bus_name, bus.d_va_deg, bus.ref_va_deg, bus.calc_va_deg)
    end
    if row.slack !== nothing
      @printf("    Slack %s dPgen=% .6g MW dQgen=% .6g MVar\n", row.slack.bus_name, row.slack.d_p_gen_MW, row.slack.d_q_gen_MVar)
    end
  end
  println()
  return nothing
end

function _scenario_file_prefix(result::ModelRunResult)::String
  scenario_tag = replace(result.scenario_name, r"[^A-Za-z0-9_.-]" => "_")
  return lowercase(result.model_name) * "_" * scenario_tag
end

function _write_summary(path::AbstractString, opt::CliOptions, ref_scenarios::Vector{For002Scenario}, branch_labels::Vector{String}, base_summaries::Vector{NamedTuple}, contingency_summaries::Vector{NamedTuple}, model_compare_file::Union{Nothing,String}, branch_b_diagnostic_csv::AbstractString, branch_b_sweep_csv::Union{Nothing,String}, worst_offenders_csv::Union{Nothing,String}, pv_bus_diagnostics_csv::Union{Nothing,String}, pv_bus_diagnostics_note::Union{Nothing,String}, active_flow_diagnostics, worst_branch_diagnostics)
  open(path, "w") do io
    println(io, "# FOR002 validation summary")
    println(io)
    println(io, "Timestamp: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
    println(io)
    println(io, "## Inputs")
    println(io, "- FOR002 reference: `", opt.for002_path, "`")
    println(io, "- MATPOWER case: `", opt.matpower_path, "`")
    println(io, "- Sparlectra builder: `", opt.builder_path, "`")
    println(io, "- Builder function: `", opt.builder_function, "`")
    println(io, "- Parsed FOR002 scenarios with data: ", length(ref_scenarios))
    println(io, "- Parsed branch labels: ", length(branch_labels))
    println(io, "- Branch b scale: ", opt.branch_b_scale)
    println(io, "- Line b scale: ", opt.line_b_scale)
    println(io, "- Transformer b scale: ", opt.trafo_b_scale)
    println(io, "- MATPOWER ratio: ", opt.matpower_ratio)
    println(io, "- Infer PV buses mode: ", opt.infer_pv_buses)
    println(io, "- Diagnose active flow: ", opt.diagnose_active_flow)
    println(io, "- Diagnose worst branches: ", opt.diagnose_worst_branches)
    println(io, "- Sweep worst branch angle: ", opt.sweep_worst_branch_angle)
    println(io, "- MATPOWER PQ generator controllers: ", opt.matpower_pq_gen_controllers)
    println(io, "- Builder PQ generator controllers: ", opt.builder_pq_gen_controllers)
    println(io, "- Builder controller handling: native-builder control is applied by clearing P(U)/Q(U) controllers on non-regulating generators after the builder returns; if none are present, the option is a no-op for that builder.")
    println(io, "- Branch charging diagnostic CSV: `", branch_b_diagnostic_csv, "`")
    if branch_b_sweep_csv !== nothing
      println(io, "- Branch charging sweep summary CSV: `", branch_b_sweep_csv, "`")
      println(io, "- Branch charging sweep objective: `max_abs_d_v_kV + 10.0 * max_abs_d_va_deg + 0.05 * max_abs_d_q_gen_MVar + 0.05 * max_abs_branch_d_q_MVar`")
    end
    if worst_offenders_csv !== nothing
      println(io, "- Worst-offender deviations CSV: `", worst_offenders_csv, "`")
    end
    if pv_bus_diagnostics_csv !== nothing
      println(io, "- PV bus diagnostics CSV: `", pv_bus_diagnostics_csv, "`")
    end
    if pv_bus_diagnostics_note !== nothing
      println(io, "- PV bus diagnostics note: ", pv_bus_diagnostics_note)
    end
    if active_flow_diagnostics !== nothing
      println(io, "- Active power balance diagnostics CSV: `", active_flow_diagnostics.active_csv, "`")
      println(io, "- Branch parameter diagnostics CSV: `", active_flow_diagnostics.branch_csv, "`")
      println(io, "- Transformer tap diagnostics CSV: `", active_flow_diagnostics.transformer_csv, "`")
    end
    if worst_branch_diagnostics !== nothing
      println(io, "- Worst branch-flow diagnostics CSV: `", worst_branch_diagnostics.flow_csv, "`")
      if worst_branch_diagnostics.sweep_csv !== nothing
        println(io, "- Worst branch angle sweep CSV: `", worst_branch_diagnostics.sweep_csv, "`")
      end
    end
    println(io)
    println(io, "## Solver settings")
    println(io, "- maxiter: ", opt.maxiter)
    println(io, "- solver_tol: ", opt.solver_tol)
    println(io, "- autodamp: ", opt.autodamp)
    println(io, "- opt_flatstart: ", opt.opt_flatstart)
    println(io, "- qlimits_enabled: ", opt.qlimits_enabled)
    println(io)
    println(io, "## Tolerances")
    println(io, "- voltage: Â±", opt.tol_v_kV, " kV / Â±", opt.tol_vm_pu, " pu")
    println(io, "- angle: Â±", opt.tol_va_deg, " deg")
    println(io, "- bus and total P/Q: Â±", opt.tol_p_MW, " MW / Â±", opt.tol_q_MVar, " MVar")
    println(io, "- branch P/Q: Â±", opt.tol_branch_p_MW, " MW / Â±", opt.tol_branch_q_MVar, " MVar")
    println(io)
    println(io, "## Base-case validation")
    println(io, "| Model | Solver status | Iterations | Residual inf | Validation |")
    println(io, "|---|---:|---:|---:|---|")
    for s in base_summaries
      println(io, "| ", s.model_name, " | ", s.status, " | ", s.iterations, " | ", something(s.residual_inf, "n/a"), " | ", s.passed ? "PASS" : "CHECK", " |")
    end
    println(io)
    for s in base_summaries
      println(io, "### ", s.model_name)
      println(io, "- Bus comparison CSV: `", s.bus_csv, "`")
      println(io, "- Branch comparison CSV: `", s.branch_csv, "`")
      println(io, "- Max |dV|: ", s.bus.max_abs_d_v_kV, " kV / ", s.bus.max_abs_d_vm_pu, " pu")
      println(io, "- Max |dVa|: ", s.bus.max_abs_d_va_deg, " deg")
      println(io, "- Max bus |dPgen| / |dQgen|: ", s.bus.max_abs_d_p_gen_MW, " MW / ", s.bus.max_abs_d_q_gen_MVar, " MVar")
      println(io, "- Max branch |dP| / |dQ|: ", s.branch.max_abs_d_p_MW, " MW / ", s.branch.max_abs_d_q_MVar, " MVar")
      println(io, "- Total loss deviation P/Q: ", s.totals.d_p_loss_MW, " MW / ", s.totals.d_q_loss_MVar, " MVar")
      println(io)
    end
    if model_compare_file !== nothing
      println(io, "## MATPOWER vs native Sparlectra builder")
      println(io, "- Bus comparison CSV: `", model_compare_file, "`")
      println(io)
    end
    if !isempty(contingency_summaries)
      println(io, "## Contingencies")
      println(io, "| Model | Scenario | Branch index out | Solver status | Iterations | dP loss MW | dQ loss MVar | max |dV| kV | max |dVa| deg |")
      println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
      for s in contingency_summaries
        println(
          io,
          "| ",
          s.model_name,
          " | ",
          s.scenario_name,
          " | ",
          something(s.branch_index_out, "n/a"),
          " | ",
          s.status,
          " | ",
          s.iterations,
          " | ",
          something(s.totals.d_p_loss_MW, "n/a"),
          " | ",
          something(s.totals.d_q_loss_MVar, "n/a"),
          " | ",
          s.bus.max_abs_d_v_kV,
          " | ",
          s.bus.max_abs_d_va_deg,
          " |",
        )
      end
      println(io)
    end
    println(io, "## Interpretation")
    println(
      io,
      "FOR002 is rounded to printed engineering values. Small residual deviations below the configured tolerances are expected. Larger systematic deviations usually indicate a conversion issue: shunt sign/model, transformer ratio side, branch status mapping, slack/generator interpretation, or branch charging treatment.",
    )
  end
  return path
end

function _load_builder_function(builder_path::AbstractString, builder_function::AbstractString, opt::CliOptions, pv_bus_diagnostic_rows::Vector{NamedTuple})
  isfile(builder_path) || throw(ArgumentError("Native Sparlectra builder file not found: $builder_path"))

  # Julia 1.12 + Revise can otherwise see the freshly included builder method
  # from an older world. Load and retrieve the function in the latest world,
  # and also call it through invokelatest when the model builder is executed.
  Base.invokelatest(include, builder_path)
  sym = Symbol(builder_function)
  fn = Base.invokelatest(() -> begin
    isdefined(Main, sym) || throw(ArgumentError("Builder function `$builder_function` was not defined by $builder_path"))
    getfield(Main, sym)
  end)
  fn isa Function || throw(ArgumentError("`$builder_function` is defined but is not a function."))
  return () -> begin
    net = Base.invokelatest(fn)
    if !opt.builder_pq_gen_controllers
      cleared = _clear_pq_gen_controllers!(net)
      cleared == 0 && @info "Builder PQ generator controller diagnostic requested, but no native-builder PQ generator P(U)/Q(U) controllers were present to clear."
    end
    append!(pv_bus_diagnostic_rows, _pv_inference_rows_and_apply!(net, opt, "builder"; apply = true))
    return net
  end
end

function _load_matpower_builder(matpower_path::AbstractString, opt::CliOptions, pv_bus_diagnostic_rows::Vector{NamedTuple})
  isfile(matpower_path) || throw(ArgumentError("MATPOWER file not found: $matpower_path"))
  mpc = Sparlectra.MatpowerIO.read_case(matpower_path; legacy_compat = true)
  return () -> begin
    net = Sparlectra.createNetFromMatPowerCase(
      mpc = mpc,
      log = false,
      flatstart = false,
      apply_bus_names = true,
      apply_branch_names = true,
      apply_branch_kind = true,
      import_for001_contingencies = true,
      matpower_ratio = opt.matpower_ratio,
      enable_pq_gen_controllers = opt.matpower_pq_gen_controllers,
    )
    append!(pv_bus_diagnostic_rows, _pv_inference_rows_and_apply!(net, opt, "matpower"; mpc = mpc, apply = true))
    return net
  end
end

function main(args = ARGS)
  opt = _parse_cli_args(collect(args))
  mkpath(opt.output_dir)

  println("FOR002 validation example")
  println("  Sparlectra version      = ", Sparlectra.version())
  println("  FOR002 reference        = ", opt.for002_path)
  println("  MATPOWER case           = ", opt.matpower_path)
  println("  Sparlectra builder      = ", opt.builder_path)
  println("  output directory        = ", opt.output_dir)
  println("  contingencies           = ", opt.run_contingencies)
  println("  branch b scale          = ", opt.branch_b_scale)
  println("  line b scale            = ", opt.line_b_scale)
  println("  transformer b scale     = ", opt.trafo_b_scale)
  println("  branch b sweep          = ", opt.sweep_branch_b)
  println("  infer PV buses          = ", opt.infer_pv_buses)
  println("  diagnose active flow    = ", opt.diagnose_active_flow)
  println("  diagnose worst branches = ", opt.diagnose_worst_branches)
  println("  sweep worst br. angle   = ", opt.sweep_worst_branch_angle)
  println("  MATPOWER ratio          = ", opt.matpower_ratio)
  println("  MATPOWER PQ controllers = ", opt.matpower_pq_gen_controllers)
  println("  builder PQ controllers  = ", opt.builder_pq_gen_controllers)
  println("  builder controller note = generic builder control clears P(U)/Q(U) controllers on non-regulating generators after build when disabled; if none exist, the option is a no-op.")
  println()

  ref_scenarios = parse_for002(opt.for002_path)
  isempty(ref_scenarios) && throw(ArgumentError("No FOR002 scenarios were parsed."))
  ref_base = ref_scenarios[1]
  branch_labels = _parse_matpower_string_vector(opt.matpower_path, "branch_name")

  pv_bus_diagnostic_rows = NamedTuple[]
  model_builders = Tuple{String,Function}[]
  if opt.run_matpower
    push!(model_builders, ("matpower", _load_matpower_builder(opt.matpower_path, opt, pv_bus_diagnostic_rows)))
  end
  if opt.run_builder
    push!(model_builders, ("builder", _load_builder_function(opt.builder_path, opt.builder_function, opt, pv_bus_diagnostic_rows)))
  end

  branch_b_diagnostic_rows = _branch_b_diagnostic_rows(opt, model_builders, branch_labels)
  branch_b_diagnostic_csv = _write_branch_b_diagnostic_csv(joinpath(opt.output_dir, "branch_b_diagnostics.csv"), branch_b_diagnostic_rows)
  _print_branch_b_diagnostic_table(branch_b_diagnostic_rows)
  println("Branch charging diagnostic CSV: ", branch_b_diagnostic_csv)
  first_pv_bus_diagnostic_rows = copy(pv_bus_diagnostic_rows)
  pv_bus_diagnostics_csv = _write_pv_bus_diagnostics_csv(joinpath(opt.output_dir, "pv_bus_diagnostics.csv"), first_pv_bus_diagnostic_rows)
  pv_bus_diagnostics_note = _pv_bus_diagnostics_note(opt, first_pv_bus_diagnostic_rows)
  println("PV bus diagnostics CSV: ", pv_bus_diagnostics_csv)
  pv_bus_diagnostics_note !== nothing && println("PV bus diagnostics note: ", pv_bus_diagnostics_note)
  println()

  branch_b_sweep_csv = nothing
  if opt.sweep_branch_b
    println("[sweep] running compact branch charging scale sweep")
    sweep_rows = _run_branch_b_sweep(opt, ref_base, model_builders, branch_labels)
    branch_b_sweep_csv = _write_branch_b_sweep_summary_csv(joinpath(opt.output_dir, "branch_b_sweep_summary.csv"), sweep_rows)
    _print_best_sweep_rows(sweep_rows)
    println("Branch charging sweep summary CSV: ", branch_b_sweep_csv)
    println()
  end

  base_results = ModelRunResult[]
  base_summaries = NamedTuple[]
  for (model_name, buildfn) in model_builders
    println("[base] running ", model_name)
    result = _run_model(model_name, buildfn, opt, "base")
    push!(base_results, result)
    prefix = _scenario_file_prefix(result)
    bus_csv = joinpath(opt.output_dir, prefix * "_bus_comparison.csv")
    branch_csv = joinpath(opt.output_dir, prefix * "_branch_comparison.csv")
    bus_cmp = _compare_bus_rows(ref_base, result, opt, bus_csv)
    branch_cmp = _compare_branch_rows(ref_base, result, branch_labels, opt, branch_csv)
    totals = _compare_totals(ref_base, result)
    passed = _base_pass(bus_cmp, branch_cmp, totals, result, opt)
    push!(base_summaries, (model_name = model_name, status = result.status, iterations = result.iterations, residual_inf = result.residual_inf, bus = bus_cmp, branch = branch_cmp, totals = totals, bus_csv = bus_csv, branch_csv = branch_csv, passed = passed))
    println("  status=", result.status, " iterations=", result.iterations, " residual=", result.residual_inf, " validation=", passed ? "PASS" : "CHECK")
  end

  model_compare_file = nothing
  if length(base_results) >= 2
    model_compare_file = joinpath(opt.output_dir, "matpower_vs_builder_bus_comparison.csv")
    _compare_two_model_runs(base_results[1], base_results[2], model_compare_file)
  end

  worst_offenders_csv = isempty(base_results) ? nothing : _write_worst_offenders_csv(joinpath(opt.output_dir, "for002_worst_offenders.csv"), ref_base, base_results, branch_labels)
  active_flow_diagnostics = opt.diagnose_active_flow ? _write_active_flow_diagnostics(opt, ref_base, base_results, branch_labels) : nothing
  _print_active_flow_console_summary(active_flow_diagnostics)
  worst_branch_diagnostics = opt.diagnose_worst_branches ? _write_worst_branch_diagnostics(opt, ref_base, base_results, model_builders, branch_labels) : nothing
  _print_worst_branch_console_summary(worst_branch_diagnostics)

  contingency_summaries = NamedTuple[]
  if opt.run_contingencies
    isempty(branch_labels) && @warn "No mpc.branch_name metadata found; outage branch mapping will be limited."
    # Use the base MATPOWER/builder branch order from the converted files.
    for ref in ref_scenarios[2:end]
      isempty(ref.buses) && continue
      branch_idx = _find_branch_index_for_outage(ref, branch_labels)
      if branch_idx === nothing
        @warn "Could not map FOR002 outage to a unique branch index; skipping scenario" scenario = ref.name outage_from = ref.outage_from outage_to = ref.outage_to outage_letter = ref.outage_letter
        continue
      end
      for (model_name, buildfn) in model_builders
        println("[contingency] running ", model_name, " / ", ref.name, " / branch ", branch_idx)
        result = _run_model(model_name, buildfn, opt, ref.name; branch_index_out = branch_idx)
        prefix = _scenario_file_prefix(result)
        bus_csv = joinpath(opt.output_dir, prefix * "_bus_comparison.csv")
        branch_csv = joinpath(opt.output_dir, prefix * "_branch_comparison.csv")
        bus_cmp = _compare_bus_rows(ref, result, opt, bus_csv)
        branch_cmp = _compare_branch_rows(ref, result, branch_labels, opt, branch_csv)
        totals = _compare_totals(ref, result)
        push!(
          contingency_summaries,
          (model_name = model_name, scenario_name = ref.name, branch_index_out = branch_idx, status = result.status, iterations = result.iterations, residual_inf = result.residual_inf, bus = bus_cmp, branch = branch_cmp, totals = totals, bus_csv = bus_csv, branch_csv = branch_csv),
        )
      end
    end
  end

  summary_path = joinpath(opt.output_dir, "for002_validation_summary.md")
  _write_summary(summary_path, opt, ref_scenarios, branch_labels, base_summaries, contingency_summaries, model_compare_file, branch_b_diagnostic_csv, branch_b_sweep_csv, worst_offenders_csv, pv_bus_diagnostics_csv, pv_bus_diagnostics_note, active_flow_diagnostics, worst_branch_diagnostics)

  println()
  println("Output files written to:")
  println("  summary                = ", summary_path)
  println("  branch b diagnostics   = ", branch_b_diagnostic_csv)
  println("  PV bus diagnostics     = ", pv_bus_diagnostics_csv)
  if active_flow_diagnostics !== nothing
    println("  active power balance  = ", active_flow_diagnostics.active_csv)
    println("  branch parameters     = ", active_flow_diagnostics.branch_csv)
    println("  transformer taps      = ", active_flow_diagnostics.transformer_csv)
  end
  if worst_branch_diagnostics !== nothing
    println("  worst branch flows    = ", worst_branch_diagnostics.flow_csv)
    worst_branch_diagnostics.sweep_csv !== nothing && println("  worst branch angle    = ", worst_branch_diagnostics.sweep_csv)
  end
  branch_b_sweep_csv !== nothing && println("  branch b sweep summary = ", branch_b_sweep_csv)
  worst_offenders_csv !== nothing && println("  worst offenders        = ", worst_offenders_csv)
  for s in base_summaries
    println("  ", s.model_name, " bus comparison     = ", s.bus_csv)
    println("  ", s.model_name, " branch comparison  = ", s.branch_csv)
  end
  model_compare_file !== nothing && println("  model comparison       = ", model_compare_file)
  if !isempty(contingency_summaries)
    println("  contingency CSV files  = ", opt.output_dir)
  end
  println("Done.")

  return (summary = summary_path, output_dir = opt.output_dir)
end

if get(ENV, "SPARLECTRA_FOR002_VALIDATION_NO_MAIN", "0") != "1"
  Base.invokelatest(getfield(@__MODULE__, :main), ARGS)
end
