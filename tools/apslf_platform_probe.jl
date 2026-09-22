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
#
# file: tools/apslf_platform_probe.jl
# purpose: prints the pieces inside the AnalyticLoadFlow solve of a 3-bus
#          ring (coefficients, Taylor sum, Pade value, the two mismatch
#          paths) so a run on another platform names the first piece that
#          differs. Run: julia --startup-file=no --project=. tools/apslf_platform_probe.jl

using Sparlectra
using LinearAlgebra
using SparseArrays
import AnalyticLoadFlow

function ring3()
    net = Net(name="apslf_ring3", baseMVA=100.0)
    addBus!(net=net, busName="B1", vn_kV=110.0, vm_pu=1.02, va_deg=0.0)
    addBus!(net=net, busName="B2", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
    addBus!(net=net, busName="B3", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
    addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B3", toBus="B1", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
    addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
    addProsumer!(net=net, busName="B2", type="GENERATOR", p=20.0, q=5.0)
    addProsumer!(net=net, busName="B3", type="LOAD", p=30.0, q=10.0)
    validate!(net=net)
    return net
end

function mismatch_paths(Y, V, spec)
    # four ways to the same injections: operator, mul! on the given Y, mul! on
    # a dense copy, and a plain loop; all four must agree
    S1 = V .* conj.(Y * V)
    I2 = zeros(ComplexF64, length(V)); mul!(I2, Y, V); S2 = V .* conj.(I2)
    I3 = zeros(ComplexF64, length(V)); mul!(I3, Matrix(Y), V); S3 = V .* conj.(I3)
    Yd = Matrix(Y)
    S4 = ComplexF64[V[i] * conj(sum(Yd[i, j] * V[j] for j in eachindex(V))) for i in eachindex(V)]
    P = i -> abs(real(S1[i]) - spec.Pspec[i])
    return (operator = maximum(P(i) for i in 2:length(V)),
            mul_given = maximum(abs.(S2 .- S1)), mul_dense = maximum(abs.(S3 .- S1)), loop = maximum(abs.(S4 .- S1)))
end

function main()
    println("Julia ", VERSION, " ", Sys.KERNEL, " ", Sys.MACHINE, " cpu ", Sys.CPU_NAME, " threads ", Threads.nthreads(), " BLAS threads ", BLAS.get_num_threads(), " ", BLAS.get_config())
    println("AnalyticLoadFlow ", pkgversion(AnalyticLoadFlow))
    model = Sparlectra.buildPfModel(ring3(); include_limits=false)
    spec = Sparlectra._apslf_spec_from_model(model)
    println("Y is ", typeof(spec.Y))
    net_nr = ring3(); runpf!(net_nr, 30, 1e-10, 0)
    vm_nr = [net_nr.nodeVec[i]._vm_pu for i in model.busIdx_net]
    println("Newton |V| ", vm_nr)
    for (label, Y) in (("sparse", spec.Y), ("dense", Matrix(spec.Y)))
        for use_pade in (true, false)
            res = AnalyticLoadFlow.solve_pf_apslf(merge(spec, (Y=Y,)); mode=:direct, order=24, use_pade=use_pade, nr_polish=false, return_coeffs=true)
            println("--- Y ", label, ", use_pade ", use_pade, ": converged ", res.converged, ", effective_mode ", res.effective_mode, ", |V| ", abs.(res.V), ", max |dVm| vs Newton ", maximum(abs.(abs.(res.V) .- vm_nr)))
            println("    mismatch paths on the returned V (P max over non-slack, then deviations of the other three paths from the operator): ", mismatch_paths(Y, res.V, spec))
            C = res.Vcoeff
            println("    coefficient magnitudes bus 2, orders 0..5: ", [abs(C[2, k]) for k in 1:6])
            taylor = [sum(C[i, :]) for i in axes(C, 1)]
            pade = [i == spec.slack ? C[i, 1] : AnalyticLoadFlow.pade_eval(collect(C[i, :]), 12, 12) for i in axes(C, 1)]
            println("    Taylor sum |V| ", abs.(taylor), ", Pade [12/12] |V| ", abs.(pade))
            println("    mismatch of the Taylor sum ", mismatch_paths(Y, taylor, spec).operator, ", of the Pade value ", mismatch_paths(Y, pade, spec).operator)
        end
    end
end

Base.invokelatest(main)
