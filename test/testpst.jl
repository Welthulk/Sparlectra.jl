using Sparlectra
using BenchmarkTools
using Logging


function test_phaseshifters()
  x_α_max = 8.0
  maxStep = 12
  neutralStep = 0
  step = 0
  δu = 0.1
  x_0 = 4.0

  myPST = SymmetricalPhaseShifter(vn=380.0, from=1, to=2, neutralStep=neutralStep, step=step, maxStep=maxStep, δu=δu, x_0=x_0, x_α_max=x_α_max)
  @debug myPST
  for i in neutralStep:maxStep
    setCurrentStep(myPST, i)
    @debug getCurrentAngle(myPST)
    @debug getCurrentX(myPST)
  end

  
  return true
end