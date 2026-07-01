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
  )
  line_b_scale_value = nothing
  trafo_b_scale_value = nothing
  run_matpower = true
  run_builder = true
  run_contingencies = false
  autodamp = true
  opt_flatstart = false
  qlimits_enabled = false

  for arg in args
    if arg == "--contingencies"
      run_contingencies = true
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

function _scenario_file_prefix(result::ModelRunResult)::String
  scenario_tag = replace(result.scenario_name, r"[^A-Za-z0-9_.-]" => "_")
  return lowercase(result.model_name) * "_" * scenario_tag
end

function _write_summary(path::AbstractString, opt::CliOptions, ref_scenarios::Vector{For002Scenario}, branch_labels::Vector{String}, base_summaries::Vector{NamedTuple}, contingency_summaries::Vector{NamedTuple}, model_compare_file::Union{Nothing,String}, branch_b_diagnostic_csv::AbstractString)
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
    println(io, "- MATPOWER PQ generator controllers: ", opt.matpower_pq_gen_controllers)
    println(io, "- Builder PQ generator controllers: ", opt.builder_pq_gen_controllers)
    println(io, "- Builder controller handling: native-builder control is applied by clearing P(U)/Q(U) controllers on non-regulating generators after the builder returns; if none are present, the option is a no-op for that builder.")
    println(io, "- Branch charging diagnostic CSV: `", branch_b_diagnostic_csv, "`")
    println(io)
    println(io, "## Solver settings")
    println(io, "- maxiter: ", opt.maxiter)
    println(io, "- solver_tol: ", opt.solver_tol)
    println(io, "- autodamp: ", opt.autodamp)
    println(io, "- opt_flatstart: ", opt.opt_flatstart)
    println(io, "- qlimits_enabled: ", opt.qlimits_enabled)
    println(io)
    println(io, "## Tolerances")
    println(io, "- voltage: ±", opt.tol_v_kV, " kV / ±", opt.tol_vm_pu, " pu")
    println(io, "- angle: ±", opt.tol_va_deg, " deg")
    println(io, "- bus and total P/Q: ±", opt.tol_p_MW, " MW / ±", opt.tol_q_MVar, " MVar")
    println(io, "- branch P/Q: ±", opt.tol_branch_p_MW, " MW / ±", opt.tol_branch_q_MVar, " MVar")
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

function _load_builder_function(builder_path::AbstractString, builder_function::AbstractString, opt::CliOptions)
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
    return net
  end
end

function _load_matpower_builder(matpower_path::AbstractString, opt::CliOptions)
  isfile(matpower_path) || throw(ArgumentError("MATPOWER file not found: $matpower_path"))
  mpc = Sparlectra.MatpowerIO.read_case(matpower_path; legacy_compat = true)
  return () -> Sparlectra.createNetFromMatPowerCase(
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
  println("  MATPOWER ratio          = ", opt.matpower_ratio)
  println("  MATPOWER PQ controllers = ", opt.matpower_pq_gen_controllers)
  println("  builder PQ controllers  = ", opt.builder_pq_gen_controllers)
  println("  builder controller note = generic builder control clears P(U)/Q(U) controllers on non-regulating generators after build when disabled; if none exist, the option is a no-op.")
  println()

  ref_scenarios = parse_for002(opt.for002_path)
  isempty(ref_scenarios) && throw(ArgumentError("No FOR002 scenarios were parsed."))
  ref_base = ref_scenarios[1]
  branch_labels = _parse_matpower_string_vector(opt.matpower_path, "branch_name")

  model_builders = Tuple{String,Function}[]
  if opt.run_matpower
    push!(model_builders, ("matpower", _load_matpower_builder(opt.matpower_path, opt)))
  end
  if opt.run_builder
    push!(model_builders, ("builder", _load_builder_function(opt.builder_path, opt.builder_function, opt)))
  end

  branch_b_diagnostic_rows = _branch_b_diagnostic_rows(opt, model_builders, branch_labels)
  branch_b_diagnostic_csv = _write_branch_b_diagnostic_csv(joinpath(opt.output_dir, "branch_b_diagnostics.csv"), branch_b_diagnostic_rows)
  _print_branch_b_diagnostic_table(branch_b_diagnostic_rows)
  println("Branch charging diagnostic CSV: ", branch_b_diagnostic_csv)
  println()

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
  _write_summary(summary_path, opt, ref_scenarios, branch_labels, base_summaries, contingency_summaries, model_compare_file, branch_b_diagnostic_csv)

  println()
  println("Output files written to:")
  println("  summary                = ", summary_path)
  println("  branch b diagnostics   = ", branch_b_diagnostic_csv)
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
