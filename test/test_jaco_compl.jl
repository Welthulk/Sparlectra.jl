using Test
using Sparlectra
include("testgrid.jl")
pq_pv = false
net = createTest5BusNet(pq_only=pq_pv)

maxIte = 20
tol    = 1e-6

ite, erg = run_complex_nr_rectangular_for_net!(net; maxiter=maxIte, tol=tol, damp=1.0, verbose=2)

if erg == 0
    calcNetLosses!(net)
    printACPFlowResults(net, 0.0, ite, tol)
else
    @warn "Rectangular complex-state PF did not converge"
end
