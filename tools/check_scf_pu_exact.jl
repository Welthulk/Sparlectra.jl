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

# file: tools/check_scf_pu_exact.jl
# purpose: checks that a MATPOWER case written to SCF with units = pu reads
#          back with the same branch and shunt per-unit values as the direct
#          import, bit for bit. Two paths per case:
#            adapter_pu: convert_case (MatpowerAdapter) -> file -> build_net
#            export_pu:  direct net -> exportSCF(units = :pu) -> build_net
#          Compared with isequal: r_pu, x_pu, b_pu, g_pu of every branch and
#          y_pu_shunt of every shunt. Prints per case and path the fields
#          that differ and exits with status 1 if any field differs.
#          Ratings and phase-shift angles are not compared: they travel in
#          SI (see docs/src/scf.md, unit system).
#
# usage:   julia --project=. tools/check_scf_pu_exact.jl [case.m ...]
#          without arguments: the MATPOWER cases shipped in data/mpower

using Sparlectra

const BRANCH_FIELDS = (:r_pu, :x_pu, :b_pu, :g_pu)

# the direct import with the configuration a run of this case would use
function direct_net(path::AbstractString, cfg)
  mpc = Sparlectra.MatpowerIO.read_case(path; legacy_compat = true)
  net = Sparlectra.createNetFromMatPowerCase(
    mpc = mpc,
    bus_shunt_model = cfg.model.bus_shunt_model,
    matpower_shift_sign = cfg.matpower.shift_sign,
    matpower_shift_unit = cfg.matpower.shift_unit,
    matpower_ratio = cfg.matpower.ratio,
    tap_changer_model = cfg.model.tap_changer_model,
    matpower_pv_voltage_source = cfg.matpower.pv_voltage_source,
    apply_branch_kind = cfg.matpower.apply_branch_kind,
    matpower_dcline_mode = cfg.matpower.matpower_dcline_mode,
  )
  return mpc, net
end

# "branch 12 x_pu" style entries for every value that is not bitwise equal
function differing_fields(A, B)::Vector{String}
  out = String[]
  length(A.branchVec) == length(B.branchVec) || return ["branch count $(length(A.branchVec)) vs $(length(B.branchVec))"]
  length(A.shuntVec) == length(B.shuntVec) || return ["shunt count $(length(A.shuntVec)) vs $(length(B.shuntVec))"]
  for i in eachindex(A.branchVec), f in BRANCH_FIELDS
    isequal(getfield(A.branchVec[i], f), getfield(B.branchVec[i], f)) || push!(out, "branch $(i) $(f)")
  end
  for i in eachindex(A.shuntVec)
    isequal(A.shuntVec[i].y_pu_shunt, B.shuntVec[i].y_pu_shunt) || push!(out, "shunt $(i) y_pu_shunt")
  end
  return out
end

function check_case(path::AbstractString, tmp::AbstractString)::Int
  cfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, path).config
  mpc, A = direct_net(path, cfg)
  fa = joinpath(tmp, "adapter_pu.scf.json")
  Sparlectra.write_scf_json(Sparlectra.convert_case(Sparlectra.MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg)), fa)
  fe = joinpath(tmp, "export_pu.scf.json")
  Sparlectra.exportSCF(A; file = fe, units = :pu)
  bad = 0
  for (label, file) in (("adapter_pu", fa), ("export_pu", fe))
    B = Sparlectra.build_net(Sparlectra.read_scf_json(file); config = cfg)
    diffs = differing_fields(A, B)
    nvals = length(A.branchVec) * length(BRANCH_FIELDS) + length(A.shuntVec)
    if isempty(diffs)
      println(basename(path), " ", label, ": all ", nvals, " values exact")
    else
      println(basename(path), " ", label, ": ", length(diffs), " of ", nvals, " values differ: ", join(first(diffs, 10), ", "), length(diffs) > 10 ? ", ..." : "")
    end
    bad += length(diffs)
  end
  return bad
end

function main(args::Vector{String})::Int
  cases = isempty(args) ? sort(filter(f -> endswith(f, ".m"), readdir(joinpath(pkgdir(Sparlectra), "data", "mpower"); join = true))) : args
  tmp = mktempdir()
  total = sum(p -> check_case(p, tmp), cases; init = 0)
  println(total == 0 ? "OK: no differing value in $(length(cases)) case(s)" : "FAIL: $(total) differing value(s)")
  return total == 0 ? 0 : 1
end

exit(Base.invokelatest(main, ARGS))
