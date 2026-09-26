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

# file: tools/dump_ybus.jl
# purpose: dump the Y-bus of every MATPOWER case under data/mpower (or of
#          the case files given as arguments) as plain text, one line per
#          nonzero entry with full-precision real and imaginary part, so
#          two checkouts can be compared entry by entry:
#            julia --project=. tools/dump_ybus.jl <outdir> [case files...]
#            diff -r <outdir_before> <outdir_after>
#          A branch-model change that must not alter the MATPOWER stamp is
#          checked with this dump before and after the change.

using Sparlectra
using SparseArrays

function dump_ybus(case_path::AbstractString, out_path::AbstractString)
  cfg = load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
  ic = Sparlectra.import_case(case_path, cfg)
  Y = createYBUS(net = ic.net, sparse = true)
  rows, cols, vals = findnz(Y)
  order = sortperm(collect(zip(rows, cols)))
  open(out_path, "w") do io
    println(io, "# ", basename(case_path), " buses=", size(Y, 1), " nnz=", length(vals))
    for k in order
      println(io, rows[k], " ", cols[k], " ", repr(real(vals[k])), " ", repr(imag(vals[k])))
    end
  end
  return length(vals)
end

function main(args)
  isempty(args) && error("usage: julia --project=. tools/dump_ybus.jl <outdir> [case files...]")
  outdir = args[1]
  mkpath(outdir)
  cases = length(args) > 1 ? args[2:end] : filter(f -> endswith(lowercase(f), ".m"), readdir(joinpath(dirname(@__DIR__), "data", "mpower"); join = true))
  for case in sort(cases)
    n = dump_ybus(case, joinpath(outdir, basename(case) * ".ybus.txt"))
    println(basename(case), ": ", n, " entries")
  end
end

Base.invokelatest(main, ARGS)
