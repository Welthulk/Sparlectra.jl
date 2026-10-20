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

# file: tools/build_scf_series.jl
# purpose: build the shipped SCF time series data/qsts/sp_case14_day.qsts.zip
#          from data/scf/sp_case14.scf.json: 96 steps of 15 minutes on
#          15 June 2026 (UTC). Loads follow a daily shape by voltage level,
#          the machines other than the reference follow the load, one line
#          is out in the afternoon. Taps and the shunt stay as in the base
#          case; moving them is the work of the controllers. The check mode
#          (default) builds into a temporary directory and compares the zip
#          entries with the shipped file; --write replaces the shipped file.
#          The random numbers come from a SplitMix64 generator written out
#          here, not from Random, so the check gives the same answer on
#          every Julia version.
#
# usage: julia --project=. tools/build_scf_series.jl [--write]

using Sparlectra

const SERIES_NAME = "sp_case14_day"
const SERIES_SEED = 20261013
const SERIES_DATE = "2026-06-15"
const SERIES_STEP_MINUTES = 15
const SERIES_STEPS = 96
const BASE_CASE = normpath(joinpath(@__DIR__, "..", "data", "scf", "sp_case14.scf.json"))
const SHIPPED_SERIES = normpath(joinpath(@__DIR__, "..", "data", "qsts", "sp_case14_day.qsts.zip"))
# one leg of the double circuit at the reference bus, out from 14:00 to 15:45
const OUTAGE_BRANCH = "Ostheim_Neuwiese_1"
const OUTAGE_STEPS = 57:64

# Daily load shapes, one value per full hour (0 to 23 h, 24 h equals 0 h),
# linear in between. Household (below 30 kV): minimum 0.45 at 04:00, morning
# peak 0.75 around 07:30, 0.70 at noon, evening peak 1.00 at 19:00.
# Commercial (30 kV and above): 0.50 at night, ramp from 06:00 to 1.00 at
# 09:00, plateau to 17:00, ramp down to 0.55 at 21:00.
const HOUSEHOLD = [0.55, 0.50, 0.47, 0.46, 0.45, 0.50, 0.60, 0.75, 0.75, 0.68, 0.66, 0.68,
                   0.70, 0.70, 0.66, 0.64, 0.68, 0.78, 0.92, 1.00, 0.97, 0.88, 0.75, 0.63]
const COMMERCIAL = [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.67, 0.83, 1.00, 1.00, 1.00,
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.89, 0.78, 0.66, 0.55, 0.52, 0.50]
# per load a factor 1 + NOISE * r, r an AR(1) series with unit variance
const NOISE = 0.02
const NOISE_AR = 0.8

# SplitMix64 (Steele, Lea, Flood 2014): a fixed sequence for a fixed seed
mutable struct SplitMix64
  state::UInt64
end

function next_u64!(r::SplitMix64)::UInt64
  r.state += 0x9e3779b97f4a7c15
  z = r.state
  z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
  z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
  return z ⊻ (z >> 31)
end

# uniform in [0, 1) from the upper 53 bits
uniform!(r::SplitMix64)::Float64 = Float64(next_u64!(r) >> 11) * 2.0^-53

# one stream per component id
stream(id::Int)::SplitMix64 = SplitMix64(UInt64(SERIES_SEED) ⊻ (UInt64(id) << 32))

# AR(1) with unit variance: uniform innovations on [-sqrt(3), sqrt(3)]
function ar1_series(r::SplitMix64, n::Int)::Vector{Float64}
  e() = (2.0 * uniform!(r) - 1.0) * sqrt(3.0)
  x = zeros(n)
  x[1] = e()
  for k in 2:n
    x[k] = NOISE_AR * x[k - 1] + sqrt(1.0 - NOISE_AR^2) * e()
  end
  return x
end

# the shape value at a minute of the day
function shape_at(points::Vector{Float64}, minute::Int)::Float64
  h = minute / 60.0
  i = floor(Int, h)
  a = points[i + 1]
  b = points[mod(i + 1, 24) + 1]
  return a + (h - i) * (b - a)
end

step_minute(k::Int)::Int = (k - 1) * SERIES_STEP_MINUTES
hhmm(minute::Int)::String = string(lpad(div(minute, 60), 2, '0'), ":", lpad(rem(minute, 60), 2, '0'))

function build_series()::SCFSeries
  base = read_scf_json(BASE_CASE)
  extra = base.sparlectra.extra
  d = base.data
  u_rated = Dict(r.id => r.u_rated for r in d.node)
  loads = sort(d.sym_load; by = r -> r.id)
  # every machine but the reference follows the load; sp_case14 has no PV
  # or wind unit, both other machines are regulated 110 kV generators
  machines = sort([r for r in d.sym_gen if get(extra[string(r.id)], "reference_pri", false) !== true]; by = r -> r.id)
  outage = only([parse(Int, id) for (id, e) in extra if get(e, "name", "") == OUTAGE_BRANCH])
  noise = Dict(r.id => ar1_series(stream(r.id), SERIES_STEPS) for r in loads)
  p_base = Dict(r.id => something(r.p_specified, 0.0) for r in loads)
  p_total = sum(values(p_base))
  steps = SCFSeriesStep[]
  ops = Vector{PatchOp}[]
  for k in 1:SERIES_STEPS
    minute = step_minute(k)
    push!(steps, SCFSeriesStep(k, string(SERIES_DATE, "T", hhmm(minute), "Z"), hhmm(minute)))
    step_ops = PatchOp[]
    weighted = 0.0
    for r in loads
      points = u_rated[r.node] < 30.0e3 ? HOUSEHOLD : COMMERCIAL
      f = round(shape_at(points, minute) * (1.0 + NOISE * noise[r.id][k]); digits = 4)
      weighted += p_base[r.id] * f
      push!(step_ops, PatchOp(op = :scale, target = :load, id = r.id, factor = f))
    end
    follow = round(weighted / p_total; digits = 4)
    for r in machines
      push!(step_ops, PatchOp(op = :scale, target = :generator, id = r.id, field = :p, factor = follow))
    end
    k in OUTAGE_STEPS && push!(step_ops, PatchOp(op = :status, target = :branch, id = outage, value = 0.0))
    push!(ops, step_ops)
  end
  provenance = Dict{String,Any}(
    "base_case" => "data/scf/sp_case14.scf.json",
    "seed" => SERIES_SEED,
    "loads" => "scale per load: daily shape by voltage level (below 30 kV household, else commercial), times 1 + 0.02 r, r an AR(1) series with a = 0.8",
    "machines" => "scale p of every machine but the reference with the load-weighted mean factor of the step",
    "outage" => string(OUTAGE_BRANCH, " out in steps ", first(OUTAGE_STEPS), " to ", last(OUTAGE_STEPS)),
  )
  return SCFSeries(name = SERIES_NAME, base = base, steps = steps, ops = ops, step_minutes = SERIES_STEP_MINUTES,
                   provenance = provenance, created_by = "tools/build_scf_series.jl")
end

# the names of the entries in which a fresh build differs from `shipped`;
# empty when they match
function check_shipped(; shipped::AbstractString = SHIPPED_SERIES)::Vector{String}
  isfile(shipped) || return ["(file missing)"]
  built = write_scf_series(joinpath(mktempdir(), basename(shipped)), build_series())
  a = Sparlectra._scf_series_read_entries(built)
  b = Sparlectra._scf_series_read_entries(shipped)
  return [name for name in Sparlectra._SCF_SERIES_ENTRIES if a[name] != b[name]]
end

function main(args::Vector{String})::Int
  unknown = filter(!=("--write"), args)
  isempty(unknown) || (println(stderr, "unknown argument(s): ", join(unknown, " "), "; usage: julia --project=. tools/build_scf_series.jl [--write]"); return 2)
  rel = relpath(SHIPPED_SERIES, normpath(joinpath(@__DIR__, "..")))
  if "--write" in args
    series = build_series()
    write_scf_series(SHIPPED_SERIES, series)
    println("written: ", rel, " (", length(series.steps), " steps, ", sum(length, series.ops), " rows, ", filesize(SHIPPED_SERIES), " bytes)")
    return 0
  end
  diffs = check_shipped()
  if isempty(diffs)
    println("check: ", rel, " matches the build")
    return 0
  end
  println("check: ", rel, " differs from the build in ", join(diffs, ", "), "; rerun with --write if the change is intended")
  return 1
end

if abspath(PROGRAM_FILE) == @__FILE__
  exit(Base.invokelatest(main, ARGS))
end
