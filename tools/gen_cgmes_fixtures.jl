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

# file: tools/gen_cgmes_fixtures.jl
# purpose: build the checked-in CGMES test deliveries under
#          data/cgmes_demo/<case>/ from the shipped cases: import the
#          case, solve it with the rectangular solver, export EQ+TP+SSH+SV
#          with writeCGMESFiles under a FIXED header stamp. The output is
#          byte-reproducible (deterministic uuid5 ids, pinned `created`), so
#          rerunning the script on an unchanged case leaves git clean. No
#          zip is written; the tests pack one at run time when they need it.
#
# usage: julia --project=. tools/gen_cgmes_fixtures.jl [sp_case14 ...]

using Sparlectra
using Dates
using Printf

# One stamp for every fixture: `created` is the only time-dependent value
# in the exported profiles, so pinning it is what makes the regeneration
# byte-identical. The tests that re-export an imported fixture use the same
# value (test/test_cgmes_export.jl, _EXPORT_STAMP).
const CGMES_FIXTURE_STAMP = DateTime(2026, 1, 1, 12, 0, 0)

const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const _FIXTURE_ROOT = normpath(joinpath(@__DIR__, "..", "data", "cgmes_demo"))

# The MATPOWER import options of the workshops and of load_fixture_net in
# the test suite (test/test_runner_helpers.jl); the fixture must be the same
# network the tests compare against.
function _load_matpower(path::AbstractString)
  return createNetFromMatPowerFile(
    filename = path,
    flatstart = false,
    enable_pq_gen_controllers = true,
    bus_shunt_model = :admittance,
    matpower_shift_sign = 1.0,
    matpower_shift_unit = :deg,
    matpower_ratio = :normal,
    tap_changer_model = :ideal,
  )
end

# case name => loader. sp_case14 carries the OLTC tap controller, sp_case118
# the Q-limited units and synchronous condensers, sp_casePST the
# phase-shifting transformer and a bus link.
const CGMES_FIXTURE_CASES = [
  "sp_case14" => () -> importSCF(joinpath(_DATA_DIR, "scf", "sp_case14.scf.json")),
  "sp_case118" => () -> _load_matpower(joinpath(_DATA_DIR, "mpower", "sp_case118.m")),
  "sp_casePST" => () -> importSCF(joinpath(_DATA_DIR, "scf", "sp_casePST.scf.json")),
]

# The exporter writes ratio-tap machinery (RatioTapChanger range, neutral,
# current step) only from `PowerTransformerTaps` on a transformer winding,
# and a voltage tap controller only when every group member carries that
# machinery. The SCF and MATPOWER importers fill the branch nameplate
# instead (has_ratio_tap, tap_step, tap_min, tap_max, tap_ratio on the
# cascade grid tap = neutral / (1 + n * step)); this maps the nameplate back
# onto the winding so the OLTC controller of sp_case14 reaches the delivery
# as a TapChangerControl. Transformer branches and net.trafos pair
# positionally, the same pairing the exporter relies on.
function _winding_taps_from_nameplate!(net::Net)
  trafo_branches = [br for br in net.branchVec if occursin("_2WT_", br.comp.cName)]
  length(trafo_branches) == length(net.trafos) || return 0
  applied = 0
  for (tf, br) in zip(net.trafos, trafo_branches)
    br.has_ratio_tap && br.tap_step > 0.0 || continue
    tf.side1.taps === nothing && tf.side2.taps === nothing || continue
    base = br.ratio
    step_of(ratio) = round(Int, (base / ratio - 1.0) / br.tap_step)
    vn1 = Sparlectra.getNodeVn(net.nodeVec[br.fromBus])
    tf.side1.taps = Sparlectra.PowerTransformerTaps(
      Vn_kV = vn1,
      step = step_of(br.tap_ratio),
      lowStep = step_of(br.tap_max),
      highStep = step_of(br.tap_min),
      neutralStep = 0,
      voltageIncrement_kV = vn1 * br.tap_step,
    )
    applied += 1
  end
  return applied
end

function build_cgmes_fixture(case::AbstractString, loader)
  net = loader()
  windings = _winding_taps_from_nameplate!(net)
  # the exporter names the profile files after the net, so the fixture
  # folder and the file prefix agree regardless of what the case file calls
  # itself
  net.name = String(case)
  # Q-limits stay off in the fixture solve. A delivery carries voltage
  # targets and the SSH operating point, not the solver's active set: a
  # unit that ended Q-pinned in the source solve re-imports as a regulated
  # PV unit at its target, the importer reads minQ/maxQ as the hull of both
  # sign readings, and the switching path then decides which units pin
  # again (measured on sp_case118: 14 pinned units, 0.012 pu against the SV
  # profile even with the source limits restored). With every regulated
  # unit holding its target, the SV profile is the solution of the
  # delivery itself, and the tests reproduce it to floating-point noise.
  ite, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular, qlimits_enabled = false)
  erg == 0 || error(string(case, ": power flow did not converge (erg = ", erg, ")"))
  outdir = joinpath(_FIXTURE_ROOT, case)
  mkpath(outdir)
  # start from an empty folder so a renamed or dropped profile never lingers
  for f in readdir(outdir; join = true)
    rm(f)
  end
  notices = String[]
  files = writeCGMESFiles(net; path = outdir, created = CGMES_FIXTURE_STAMP, notices = notices)
  total = 0
  println(case, ": ", length(net.nodeVec), " buses, ", length(net.branchVec), " branches, solved in ", ite, " iterations, ", windings, " winding(s) with ratio-tap machinery from the branch nameplate")
  for f in files
    sz = filesize(f)
    total += sz
    @printf("  %-22s %9d bytes\n", basename(f), sz)
  end
  @printf("  %-22s %9d bytes\n", "total", total)
  for n in notices
    println("  exporter notice: ", n)
  end
  return total
end

function main(args)
  selected = isempty(args) ? first.(CGMES_FIXTURE_CASES) : String.(args)
  for (case, loader) in CGMES_FIXTURE_CASES
    case in selected || continue
    build_cgmes_fixture(case, loader)
  end
  unknown = setdiff(selected, first.(CGMES_FIXTURE_CASES))
  isempty(unknown) || error(string("unknown fixture case(s): ", join(unknown, ", ")))
  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(main, ARGS)
end
