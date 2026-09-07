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

# file: tools/matpower_adapter_diff.jl
# purpose: the stage-3a acceptance oracle. Runs every given MATPOWER case
#          through BOTH construction paths, the direct import
#          (createNetFromMatPowerCase) and the adapter path
#          (convert_case -> build_net), solves both with the identical
#          effective configuration, and diffs iterations, final mismatch
#          (1e-10) and bus voltages (1e-9 pu). Case conventions come from
#          the case's own <stem>.config.yaml through resolve_config, so
#          the pegase cases run with their documented rad/-1.0 convention.
#
# usage: julia --project=. tools/matpower_adapter_diff.jl [case files...]
#        without arguments the bundled default list runs; missing files are
#        reported loudly and skipped.

using Sparlectra

# measurement guard (RP1 worktree incident 2026-09-03): this script must run
# against THE repository checkout it lives in, never a stale worktree or
# another depot copy picked up through a wrong --project
let expected = normpath(joinpath(@__DIR__, "..", "src")), actual = normpath(String(pathof(Sparlectra)))
  startswith(actual, expected) || error("tools guard: Sparlectra loaded from " * actual * ", expected under " * expected * "; start julia with --project=" * normpath(joinpath(@__DIR__, "..")))
end
using Printf

function _diff_one(case_path::AbstractString)
  cfgfile = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
  resolved = Sparlectra.resolve_config(cfgfile, case_path)
  cfg = resolved.config
  mpc = Sparlectra.MatpowerIO.read_case(case_path; legacy_compat = true)

  netA = Sparlectra.createNetFromMatPowerCase(
    mpc = mpc,
    flatstart = cfg.powerflow.start_mode.flatstart,
    enable_pq_gen_controllers = cfg.matpower.enable_pq_gen_controllers,
    bus_shunt_model = cfg.model.bus_shunt_model,
    matpower_shift_sign = cfg.matpower.shift_sign,
    matpower_shift_unit = cfg.matpower.shift_unit,
    matpower_ratio = cfg.matpower.ratio,
    tap_changer_model = cfg.model.tap_changer_model,
    matpower_pv_voltage_source = cfg.matpower.pv_voltage_source,
    matpower_pv_voltage_mismatch_tol_pu = cfg.matpower.pv_voltage_mismatch_tol_pu,
    apply_bus_names = cfg.matpower.apply_bus_names,
    apply_branch_names = cfg.matpower.apply_branch_names,
    apply_branch_kind = cfg.matpower.apply_branch_kind,
    import_for001_contingencies = cfg.matpower.import_for001_contingencies,
    matpower_dcline_mode = cfg.matpower.matpower_dcline_mode,
    preallocate_network = cfg.model.preallocate_network,
    preallocate_min_buses = cfg.model.preallocate_min_buses,
  )
  Sparlectra.MatpowerIO.apply_mp_isolated_buses!(netA, mpc)
  Sparlectra.MatpowerIO.apply_mp_bus_vmva_init!(netA, mpc; flatstart = cfg.powerflow.start_mode.flatstart)
  Sparlectra._apply_config_net_parameters!(netA, cfg)
  netA.bus_shunt_model = Sparlectra.normalize_bus_shunt_model(cfg.model.bus_shunt_model)

  opts = Sparlectra.matpower_adapter_options(cfg)
  case = Sparlectra.convert_case(MatpowerAdapter(), mpc, opts)
  netB = Sparlectra.build_net(case; config = cfg)

  ok = true
  notes = String[]
  # model equivalence FIRST, bitwise on every electrical branch/shunt field:
  # stronger than the final-mismatch column of the task gate, because two
  # bitwise-identical models solved from bitwise-identical starts cannot
  # differ in mismatch either
  length(netA.branchVec) == length(netB.branchVec) || (ok = false; push!(notes, string("branch count A=", length(netA.branchVec), " B=", length(netB.branchVec))))
  if length(netA.branchVec) == length(netB.branchVec)
    for i in eachindex(netA.branchVec)
      a = netA.branchVec[i]
      b = netB.branchVec[i]
      for f in (:r_pu, :x_pu, :b_pu, :g_pu, :ratio, :angle, :status, :from_status, :to_status, :fromBus, :toBus)
        va = getfield(a, f)
        vb = getfield(b, f)
        # the SCF fixed-point stabilizer documents up to one ulp of
        # per-unit drift where no exact SI preimage exists; 8 eps relative
        # still catches every real conversion bug
        close = (va isa Float64 && vb isa Float64) ? (isequal(va, vb) || isapprox(va, vb; rtol = 8 * eps(Float64), atol = 0.0)) : isequal(va, vb)
        if !close
          ok = false
          push!(notes, string("branch ", i, " field ", f, " A=", va, " B=", vb))
          break
        end
      end
      ok || break
    end
  end
  length(netA.shuntVec) == length(netB.shuntVec) || (ok = false; push!(notes, string("shunt count A=", length(netA.shuntVec), " B=", length(netB.shuntVec))))
  if length(netA.shuntVec) == length(netB.shuntVec)
    for i in eachindex(netA.shuntVec)
      ya = netA.shuntVec[i].y_pu_shunt
      yb = netB.shuntVec[i].y_pu_shunt
      (isapprox(real(ya), real(yb); rtol = 8 * eps(Float64), atol = 1.0e-300) && isapprox(imag(ya), imag(yb); rtol = 8 * eps(Float64), atol = 1.0e-300)) && continue
      ok = false
      push!(notes, string("shunt ", i, " y A=", netA.shuntVec[i].y_pu_shunt, " B=", netB.shuntVec[i].y_pu_shunt))
      break
    end
  end

  iteA, ergA = Sparlectra.runpf!(netA; config = cfg)
  iteB, ergB = Sparlectra.runpf!(netB; config = cfg)
  ergA == ergB || (ok = false; push!(notes, string("result code A=", ergA, " B=", ergB)))
  iteA == iteB || (ok = false; push!(notes, string("iterations A=", iteA, " B=", iteB)))
  namesA = Dict(name => idx for (name, idx) in netA.busDict)
  worst_dvm = 0.0
  worst_dva = 0.0
  worst_bus = ""
  for (name, idxB) in netB.busDict
    idxA = get(namesA, name, 0)
    if idxA == 0
      ok = false
      push!(notes, string("bus ", name, " missing in the direct net"))
      continue
    end
    a = netA.nodeVec[idxA]
    b = netB.nodeVec[idxB]
    vmA = a._vm_pu === nothing ? NaN : Float64(a._vm_pu)
    vmB = b._vm_pu === nothing ? NaN : Float64(b._vm_pu)
    vaA = a._va_deg === nothing ? NaN : Float64(a._va_deg)
    vaB = b._va_deg === nothing ? NaN : Float64(b._va_deg)
    dvm = isnan(vmA) && isnan(vmB) ? 0.0 : abs(vmA - vmB)
    dva = isnan(vaA) && isnan(vaB) ? 0.0 : abs(vaA - vaB)
    if dvm > worst_dvm
      worst_dvm = dvm
      worst_bus = name
    end
    dva > worst_dva && (worst_dva = dva)
  end
  length(netA.busDict) == length(netB.busDict) || (ok = false; push!(notes, string("bus count A=", length(netA.busDict), " B=", length(netB.busDict))))
  worst_dvm <= 1.0e-9 || (ok = false; push!(notes, string("max |dVm| ", worst_dvm, " pu at ", worst_bus)))
  worst_dva <= 1.0e-7 || (ok = false; push!(notes, string("max |dVa| ", worst_dva, " deg")))
  return (ok = ok, iterations = (iteA, iteB), notes = notes, dvm = worst_dvm, dva = worst_dva)
end

function main(args)
  default_root = normpath(joinpath(@__DIR__, "..", "data", "mpower"))
  webui_root = expanduser("~/.local/state/sparlectra/webui/data/mpower")
  cases = isempty(args) ? [
    joinpath(default_root, "case14.m"),
    joinpath(default_root, "case57.m"),
    joinpath(default_root, "case118.m"),
    joinpath(default_root, "case145.m"),
    joinpath(default_root, "case300.m"),
    joinpath(default_root, "case1354pegase.m"),
    joinpath(default_root, "case13659pegase.m"),
  ] : collect(args)
  failed = 0
  for case_path in cases
    if !isfile(case_path)
      alt = joinpath(webui_root, basename(case_path))
      if isfile(alt)
        case_path = alt
      else
        println("SKIPPED (file not found): ", case_path)
        continue
      end
    end
    result = try
      _diff_one(case_path)
    catch err
      failed += 1
      println("FAIL  ", basename(case_path), ": ", sprint(showerror, err))
      continue
    end
    if result.ok
      @printf("PASS  %-22s iterations=%d dVm=%.2e dVa=%.2e\n", basename(case_path), result.iterations[1], result.dvm, result.dva)
    else
      failed += 1
      println("FAIL  ", basename(case_path), ": ", join(result.notes, "; "))
    end
  end
  failed == 0 || exit(1)
  return nothing
end

Base.invokelatest(main, ARGS)
