# Copyright 2023–2025 Udo Schmitz
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

module MatpowerIO

export MatpowerCase, read_case, read_case_m, read_case_julia, build_ybus_matpower

using LinearAlgebra

"""
Container for a MATPOWER case (case format v2/v2-ish).
All matrices are stored as Float64 matrices, names as Vector{String} when available.
"""
struct MatpowerCase
  name::String
  baseMVA::Float64
  bus::Matrix{Float64}
  gen::Matrix{Float64}
  branch::Matrix{Float64}
  gencost::Union{Nothing,Matrix{Float64}}
  bus_name::Union{Nothing,Vector{String}}
end

"""
    legacy_sort_bus(mpc::MatpowerCase) -> MatpowerCase

Mimics the legacy `casefileparser` behavior:
- Sorts `bus` rows by BUS_I (column 1) ascending.
- Does NOT reorder `gen` or `branch`.
This keeps results compatible with the old importer.
"""
function legacy_sort_bus(mpc::MatpowerCase)::MatpowerCase
  bus = mpc.bus
  perm = sortperm(bus[:, 1])  # sort by BUS_I
  bus_sorted = bus[perm, :]
  return MatpowerCase(mpc.name, mpc.baseMVA, bus_sorted, mpc.gen, mpc.branch, mpc.gencost, mpc.bus_name)
end

"""
    normalize_branch_tap!(branch)

MATPOWER semantics:
- column 9 (TAP) == 0 means "no transformer tap specified" -> treat as 1.0.

We normalize data early so downstream code sees a consistent ratio.
"""
function normalize_branch_tap!(branch::AbstractMatrix{<:Real})
  ncol = size(branch, 2)
  ncol < 9 && return branch
  anz = size(branch, 1)
  @inbounds for e = 1:anz
    if branch[e, 9] == 0.0
      branch[e, 9] = 1.0
    end
  end
  return branch
end

"""
    read_case(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

Dispatches based on file extension:
- `.m`  => `read_case_m`
- `.jl` => `read_case_julia` (expects file returns `MatpowerCase` or `NamedTuple`)
"""
function read_case(path::AbstractString; legacy_compat::Bool = true)
  lpath = lowercase(path)
  mpc = if endswith(lpath, ".m")
    read_case_m(path; legacy_compat = legacy_compat)
  elseif endswith(lpath, ".jl")
    read_case_julia(path; legacy_compat = legacy_compat)
  else
    error("Unsupported case file extension. Expected .m or .jl, got: $path")
  end
  return mpc
end

"""
    read_case_julia(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

The `.jl` file should either:
(A) return a `MatpowerCase`, OR
(B) return a `NamedTuple` with fields: `baseMVA, bus, gen, branch` and optional `gencost, bus_name, name`.
"""
function read_case_julia(path::AbstractString; legacy_compat::Bool = true)
  obj = include(path)

  mpc = if obj isa MatpowerCase
    obj
  elseif obj isa NamedTuple
    name = get(obj, :name, splitext(basename(path))[1])
    baseMVA = Float64(obj.baseMVA)
    bus = Matrix{Float64}(obj.bus)
    gen = Matrix{Float64}(obj.gen)
    branch = Matrix{Float64}(obj.branch)
    normalize_branch_tap!(branch)

    gencost = haskey(obj, :gencost) ? (obj.gencost === nothing ? nothing : Matrix{Float64}(obj.gencost)) : nothing
    bus_name = haskey(obj, :bus_name) ? (obj.bus_name === nothing ? nothing : Vector{String}(obj.bus_name)) : nothing
    MatpowerCase(String(name), baseMVA, bus, gen, branch, gencost, bus_name)
  else
    error("Julia case file must return MatpowerCase or NamedTuple, got: $(typeof(obj))")
  end

  return legacy_compat ? legacy_sort_bus(mpc) : mpc
end

"""
    read_case_m(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

Best-effort parser for typical MATPOWER case files (*.m).
Parses:
- `mpc.baseMVA`
- `mpc.bus`
- `mpc.gen`
- `mpc.branch`
- optional: `mpc.gencost`
- optional: `mpc.bus_name`
"""
function read_case_m(path::AbstractString; legacy_compat::Bool = true)
  txt = read(path, String)
  name = splitext(basename(path))[1]

  # Remove MATLAB comments (%) but keep line structure
  lines = split(txt, '\n')
  lines_nc = String[]
  for ln in lines
    i = findfirst('%', ln)
    push!(lines_nc, i === nothing ? ln : ln[1:prevind(ln, i)])
  end
  txt = join(lines_nc, "\n")

  baseMVA = parse_baseMVA(txt)

  # fixed widths to match MATPOWER v2 tables
  bus    = parse_matrix_block(txt, "mpc.bus"; ncols = 13)
  gen    = parse_matrix_block(txt, "mpc.gen"; ncols = 21)
  branch = parse_matrix_block(txt, "mpc.branch"; ncols = 13)
  normalize_branch_tap!(branch)

  gencost = try_parse_matrix_block(txt, "mpc.gencost")
  bus_name = try_parse_bus_name(txt)

  mpc = MatpowerCase(name, baseMVA, bus, gen, branch, gencost, bus_name)
  return legacy_compat ? legacy_sort_bus(mpc) : mpc
end

function parse_baseMVA(txt::String)
  m = match(r"mpc\.baseMVA\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*;", txt)
  m === nothing && error("Could not find `mpc.baseMVA = ...;` in MATPOWER file.")
  return parse(Float64, m.captures[1])
end

function try_parse_matrix_block(txt::String, key::String)
  try
    return parse_matrix_block(txt, key)
  catch
    return nothing
  end
end

function parse_matrix_block(txt::String, key::String; ncols::Int = 0)
  re = Regex(replace(key, "." => "\\.") * raw"\s*=\s*\[\s*(.*?)\s*\]\s*;", "s")
  m = match(re, txt)
  m === nothing && error("Could not find matrix block `$key = [ ... ];`")
  body = m.captures[1]

  if ncols > 0
    return parse_numeric_matrix_ncols(body, key; ncols = ncols)
  else
    return parse_numeric_matrix(body, key)
  end
end

# best-effort: auto width; IMPORTANT: fill with 0.0 (not NaN)
function parse_numeric_matrix(body::AbstractString, key::AbstractString)
  raw_rows = split(body, ';')
  rows = Vector{Vector{Float64}}()
  maxcols = 0

  for rr in raw_rows
    s = strip(rr)
    isempty(s) && continue
    s = replace(s, '\t' => ' ')
    s = replace(s, r"\s+" => " ")
    toks = split(s, ' ')

    vals = Float64[]
    for t in toks
      tt = strip(t)
      isempty(tt) && continue
      push!(vals, parse(Float64, tt))
    end

    maxcols = max(maxcols, length(vals))
    push!(rows, vals)
  end

  isempty(rows) && error("Matrix `$key` appears empty or could not be parsed.")

  M = zeros(Float64, length(rows), maxcols)
  for (i, r) in enumerate(rows)
    for (j, val) in enumerate(r)
      M[i, j] = val
    end
  end
  return M
end

# fixed width path: pad/truncate to exactly ncols with zeros
function parse_numeric_matrix_ncols(body::AbstractString, key::AbstractString; ncols::Int)
  raw_rows = split(body, ';')
  rows = Vector{Vector{Float64}}()

  for rr in raw_rows
    s = strip(rr)
    isempty(s) && continue
    s = replace(s, '\t' => ' ')
    s = replace(s, r"\s+" => " ")
    toks = split(s, ' ')

    vals = Float64[]
    for t in toks
      tt = strip(t)
      isempty(tt) && continue
      push!(vals, parse(Float64, tt))
    end

    if length(vals) < ncols
      append!(vals, zeros(Float64, ncols - length(vals)))
    elseif length(vals) > ncols
      vals = vals[1:ncols]
    end

    push!(rows, vals)
  end

  isempty(rows) && error("Matrix `$key` appears empty or could not be parsed.")

  M = zeros(Float64, length(rows), ncols)
  for (i, r) in enumerate(rows)
    @inbounds M[i, :] .= r
  end
  return M
end

function try_parse_bus_name(txt::String)
  re = r"mpc\.bus_name\s*=\s*\{\s*(.*?)\s*\}\s*;"s
  m = match(re, txt)
  m === nothing && return nothing
  body = m.captures[1]

  names = String[]
  for mm in eachmatch(r"'([^']*)'", body)
    push!(names, mm.captures[1])
  end
  isempty(names) ? nothing : names
end

"""
    build_ybus_matpower(bus, branch, baseMVA) -> Matrix{ComplexF64}

MATPOWER-style Ybus stamping (π-model + tap/shift + bus shunts):
- series r/x
- line charging b split as b/2 on each end
- off-nominal tap ratio and phase shift on from-side
- bus shunts GS/BS added on diagonal: (GS + j*BS)/baseMVA
"""
function build_ybus_matpower(bus::AbstractMatrix{<:Real}, branch::AbstractMatrix{<:Real}, baseMVA::Float64)
  nbus = size(bus, 1)
  Y = zeros(ComplexF64, nbus, nbus)

  # bus columns (MATPOWER):
  # BUS_I=1, BUS_TYPE=2, PD=3, QD=4, GS=5, BS=6, ...
  # add shunts: (GS + j*BS)/baseMVA
  for i = 1:nbus
    Gs = bus[i, 5]
    Bs = bus[i, 6]
    if Gs != 0.0 || Bs != 0.0
      Y[i, i] += (Gs + 1im*Bs) / baseMVA
    end
  end

  # branch columns (MATPOWER):
  # F_BUS=1, T_BUS=2, BR_R=3, BR_X=4, BR_B=5, RATEA=6, RATEB=7, RATEC=8,
  # TAP=9, SHIFT=10, BR_STATUS=11, ...
  nbranch = size(branch, 1)
  for e = 1:nbranch
    status = (size(branch, 2) >= 11) ? branch[e, 11] : 1.0
    status == 0 && continue

    f = Int(branch[e, 1])
    t = Int(branch[e, 2])
    r = branch[e, 3]
    x = branch[e, 4]
    b = branch[e, 5]

    z = r + 1im*x
    y = inv(z)

    # line charging b -> j*b/2 at both ends
    ysh = 1im * b / 2

    tap = (size(branch, 2) >= 9) ? branch[e, 9] : 1.0
    shift_deg = (size(branch, 2) >= 10) ? branch[e, 10] : 0.0

    # enforce MATPOWER semantics again (safe)
    if tap == 0.0
      tap = 1.0
    end

    shift = shift_deg * (pi / 180)
    a = tap * cis(shift)  # complex tap on from side

    # stamping with off-nominal tap
    Y[f, f] += (y + ysh) / (a * conj(a))
    Y[t, t] += (y + ysh)
    Y[f, t] += -y / conj(a)
    Y[t, f] += -y / a
  end

  return Y
end



end # module MatpowerIO
