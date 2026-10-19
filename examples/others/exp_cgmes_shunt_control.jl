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

# file: examples/others/exp_cgmes_shunt_control.jl
# purpose: cgmes_import.shunt_control on the ENTSO-E SmallGrid: the banks
#          with a voltage RegulatingControl become switched-bank regulators
#          that take over the imported shunts; a raised target at one bank
#          makes it switch blocks in the outer control loop.
#          julia --project=. examples/others/exp_cgmes_shunt_control.jl

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const SMALLGRID = sort!(filter(f -> endswith(f, ".zip"), readdir(joinpath(dirname(dirname(pathof(Sparlectra))), "test", "fixtures", "cgmes_smallgrid"); join = true)))

function main()
  print_example_banner("examples/others/exp_cgmes_shunt_control.jl", "CGMES bank regulators: ShuntCompensator plus RegulatingControl as switched banks")

  # 1. the import with and without the switch
  off = redirect_stdout(() -> importCGMES(path = SMALLGRID, name = "smallgrid"), devnull)
  on = redirect_stdout(() -> importCGMES(path = SMALLGRID, name = "smallgrid", shunt_control = true), devnull)
  banks = Sparlectra._shunt_controllers(on.net)
  println("shunts ", length(on.net.shuntVec), ", bank regulators without the switch ", length(Sparlectra._shunt_controllers(off.net)), ", with it ", length(banks))
  foreach(m -> occursin("stays a fixed shunt", m) && println("  ", m), on.messages)

  # 2. the run with the controllers: the delivery is solved within every band
  redirect_stdout(() -> run_sparlectra(net = on.net, config = SparlectraConfig()), devnull)
  println("control status ", latest_control_result(on.net).status)

  # 3. one bank with a target 0.02 pu higher switches blocks in
  net = redirect_stdout(() -> importCGMES(path = SMALLGRID, name = "smallgrid", shunt_control = true), devnull).net
  bank = first(Sparlectra._shunt_controllers(net))
  start = bank.bs_mvar
  bank.target_vm_pu += 0.02
  redirect_stdout(() -> run_sparlectra(net = net, config = SparlectraConfig()), devnull)
  @printf("%s at %s: target %.4f pu, %.3f MVAr -> %.3f MVAr (%d blocks of %.3f MVAr), voltage %.4f pu, status %s\n", bank.name, bank.bus, bank.target_vm_pu, start, bank.bs_mvar,
    round(Int, bank.bs_mvar / bank.step_mvar), bank.step_mvar, get_bus_vm_pu(net, bank.bus), latest_control_result(net).status)
  return (banks = length(banks), moved = bank.bs_mvar != start)
end

result = run_example(main)
println()
println("bank regulators: ", result.banks, ", the raised target moved its bank: ", result.moved)
