using Sparlectra
using BenchmarkTools
using Logging


function test_phaseshifters()
  myPST = SymmetricalPhaseShifter(vn=380.0, from=1, to=2, neutralStep=0, step=6, maxStep=12, δu=0.1, x_0=4.0, x_α_max=8.0)
  @show myPST
  return true
end