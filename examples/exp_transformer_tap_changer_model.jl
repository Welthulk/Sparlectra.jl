# Copyright 2023–2026 Udo Schmitz
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

using Sparlectra

"""
    main()

Import the same minimal MATPOWER transformer branch once with the default
`:ideal` tap-changer model and once with `:impedance_correction`, and compare
the resulting series impedance and power-flow solution.

`:ideal` keeps the tap changer free of impedance feedback: R/X stay at the
imported neutral-position values regardless of the off-nominal tap ratio.
`:impedance_correction` re-refers R/X through the tapped winding, scaling
them with the squared magnitude of the regulating vector `|1 + f·e^(jφ)|²`.
Both models are implemented centrally in `calcTapCorrectedRX` /
`calcTapImpedanceCorrectionFactor` (`src/equicircuit.jl`) and are read by both
the MATPOWER and the native DTF importer through
`transformer.tap_changer_model`.
"""
function main()
  mpc = (
    name = "tap_changer_model_demo",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      2 1 20.0 10.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
    ],
    gen = [1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0;],
    # branch: [fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax]
    # ratio = 0.93 models a transformer running 7% off the neutral tap position.
    branch = [1 2 0.01 0.10 0.0 9999.0 0.0 0.0 0.93 0.0 1 -60.0 60.0;],
    gencost = nothing,
    bus_name = nothing,
  )

  net_ideal = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true, tap_changer_model = :ideal)
  net_corrected = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true, tap_changer_model = :impedance_correction)

  ite_ideal, status_ideal = runpf!(net_ideal, 30, 1e-8, 0)
  ite_corrected, status_corrected = runpf!(net_corrected, 30, 1e-8, 0)

  branch_ideal = only(net_ideal.branchVec)
  branch_corrected = only(net_corrected.branchVec)

  return (
    ratio = 0.93,
    r_pu_ideal = branch_ideal.r_pu,
    x_pu_ideal = branch_ideal.x_pu,
    r_pu_impedance_correction = branch_corrected.r_pu,
    x_pu_impedance_correction = branch_corrected.x_pu,
    correction_factor = branch_corrected.r_pu / branch_ideal.r_pu,
    converged_ideal = status_ideal == 0,
    converged_impedance_correction = status_corrected == 0,
    iterations_ideal = ite_ideal,
    iterations_impedance_correction = ite_corrected,
    vm_bus2_ideal = getNodeVm(net_ideal.nodeVec[2]),
    vm_bus2_impedance_correction = getNodeVm(net_corrected.nodeVec[2]),
  )
end

if abspath(PROGRAM_FILE) == @__FILE__
  result = Base.invokelatest(main)
  println("Off-nominal tap ratio: ", result.ratio)
  println("ideal:                R=", result.r_pu_ideal, " pu  X=", result.x_pu_ideal, " pu  converged=", result.converged_ideal, " (", result.iterations_ideal, " it.)  Vm(bus2)=", result.vm_bus2_ideal, " pu")
  println("impedance_correction: R=", result.r_pu_impedance_correction, " pu  X=", result.x_pu_impedance_correction, " pu  converged=", result.converged_impedance_correction, " (", result.iterations_impedance_correction, " it.)  Vm(bus2)=", result.vm_bus2_impedance_correction, " pu")
  println("R/X correction factor: ", result.correction_factor)
end
