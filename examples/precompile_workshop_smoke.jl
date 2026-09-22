# file: precompile_workshop_smoke.jl
# purpose: Standalone-Version des Workload-Blocks zum Ausprobieren in REPL
#          oder als Skript. Gleicher Pfad wie die Workshop-Notebooks:
#          Netz programmatisch bauen, run_sparlectra mit rechteckigem und
#          APSLF-Solver, Statusabfrage, APSLF-Start vor NR.
#          Außerhalb des Moduls sind ApslfConfig, ApslfStartConfig und
#          rectangular_pf_status nicht exportiert, daher qualifiziert.
#
# Aufruf: julia --project=. precompile_workshop_smoke.jl
# Die @time-Zeilen zeigen die Kompilierzeit des ersten Aufrufs; nach dem
# Einbau des Blocks in src/build/precompile.jl sollten sie deutlich fallen.

using Sparlectra

function build_ring3()
    net = Net(name="precompile_ring3", baseMVA=100.0)
    addBus!(net=net, busName="B1", vn_kV=110.0, vm_pu=1.02, va_deg=0.0)
    addBus!(net=net, busName="B2", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
    addBus!(net=net, busName="B3", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
    addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B3", toBus="B1", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
    addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
    addProsumer!(net=net, busName="B2", type="GENERATOR", p=20.0, q=5.0)
    addProsumer!(net=net, busName="B3", type="LOAD", p=30.0, q=10.0)
    ok, msg = validate!(net=net)
    ok || error("Network validation failed: $msg")
    return net
end

quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
cfg_nr = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false), output=quiet)
cfg_ap = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=12)), output=quiet)
cfg_hyb = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, apslf_start=Sparlectra.ApslfStartConfig(enabled=true, order=12)), output=quiet)

println("first call (includes compilation):")
@time r_nr = run_sparlectra(net=build_ring3(), config=cfg_nr)
@time r_ap = run_sparlectra(net=build_ring3(), config=cfg_ap)
@time r_hyb = run_sparlectra(net=build_ring3(), config=cfg_hyb)

println("second call (compiled):")
@time run_sparlectra(net=build_ring3(), config=cfg_nr)
@time run_sparlectra(net=build_ring3(), config=cfg_ap)
@time run_sparlectra(net=build_ring3(), config=cfg_hyb)

st_nr = Sparlectra.rectangular_pf_status(r_nr.net)
st_ap = Sparlectra.rectangular_pf_status(r_ap.net)
println("NR    : ", r_nr.outcome, ", ", r_nr.iterations, " iterations, mismatch ", r_nr.final_mismatch)
println("APSLF : ", r_ap.outcome, ", ", r_ap.iterations, " pass(es),   mismatch ", r_ap.final_mismatch)
println("APSLF : ", st_ap.apslf_convergence_line)
println("hybrid: ", r_hyb.outcome, ", ", r_hyb.iterations, " Newton iterations")
println("Vm    : ", [n._vm_pu for n in r_ap.net.nodeVec])