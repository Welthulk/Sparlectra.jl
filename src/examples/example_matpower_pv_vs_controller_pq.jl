# Example: MATPOWER generator handling as PV vs constant P(U)/Q(U)-controlled PQ.
#
# Goal:
# - Case A: keep MATPOWER BUS_TYPE for generator buses (classic PV behavior)
# - Case B: force non-slack generator buses to PQ (BUS_TYPE=1), which maps
#           MATPOWER generator limits to constant P(U)/Q(U) controllers.
#
# This demonstrates how the importer treats generators that are not specified
# as PV buses.

using Sparlectra
using Printf
using Logging

const MPOWER_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "mpower"))

function ensure_case_available!(case::AbstractString; overwrite::Bool = false)
  mkpath(MPOWER_DIR)
  local_m = joinpath(MPOWER_DIR, case)
  if !isfile(local_m) || overwrite
    try
      url = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/$(case)"
      Sparlectra.FetchMatpowerCase.ensure_matpower_case(url = url, outdir = MPOWER_DIR, to_jl = false, overwrite = overwrite)
    catch err
      @warn "Could not download $(case). Falling back to inline case9 data." error = err
      return nothing
    end
  end
  return local_m
end

function case9_inline()
  bus = Float64[
    1 3 0 0 0 0 1 1 0 345 1 1.1 0.9
    2 2 0 0 0 0 1 1 0 345 1 1.1 0.9
    3 2 0 0 0 0 1 1 0 345 1 1.1 0.9
    4 1 0 0 0 0 1 1 0 345 1 1.1 0.9
    5 1 90 30 0 0 1 1 0 345 1 1.1 0.9
    6 1 0 0 0 0 1 1 0 345 1 1.1 0.9
    7 1 100 35 0 0 1 1 0 345 1 1.1 0.9
    8 1 0 0 0 0 1 1 0 345 1 1.1 0.9
    9 1 125 50 0 0 1 1 0 345 1 1.1 0.9
  ]
  gen = Float64[
    1 72.3 27.03 300 -300 1 100 1 250 10 0 0 0 0 0 0 0 0 0 0 0
    2 163.0 6.54 300 -300 1 100 1 300 10 0 0 0 0 0 0 0 0 0 0 0
    3 85.0 -10.95 300 -300 1 100 1 270 10 0 0 0 0 0 0 0 0 0 0 0
  ]
  branch = Float64[
    1 4 0 0.0576 0 250 250 250 0 0 1 -360 360
    4 5 0.017 0.092 0.158 250 250 250 0 0 1 -360 360
    5 6 0.039 0.17 0.358 150 150 150 0 0 1 -360 360
    3 6 0 0.0586 0 300 300 300 0 0 1 -360 360
    6 7 0.0119 0.1008 0.209 150 150 150 0 0 1 -360 360
    7 8 0.0085 0.072 0.149 250 250 250 0 0 1 -360 360
    8 2 0 0.0625 0 250 250 250 0 0 1 -360 360
    8 9 0.032 0.161 0.306 250 250 250 0 0 1 -360 360
    9 4 0.01 0.085 0.176 250 250 250 0 0 1 -360 360
  ]
  return (name = "case9_inline", baseMVA = 100.0, bus = bus, gen = gen, branch = branch)
end

function generator_summary(net::Net)
  println("Generator setup:")
  println(" bus | regulated(PV) | has P(U) | has Q(U) | P [MW] | Q [MVAr] | Vm set")
  println("-----+---------------+----------+----------+--------+----------+-------")
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bidx = getPosumerBusIndex(ps)
    @printf(" %3d | %13s | %8s | %8s | %6.2f | %8.2f | %5s\n", bidx, string(ps.isRegulated), string(has_pu_controller(ps)), string(has_qu_controller(ps)), isnothing(ps.pVal) ? 0.0 : ps.pVal, isnothing(ps.qVal) ? 0.0 : ps.qVal, isnothing(ps.vm_pu) ? "-" : @sprintf("%.3f", ps.vm_pu),)
  end
  println("")
end

function force_generator_buses_to_pq!(mpc)
  gen_buses = unique(Int(row[1]) for row in eachrow(mpc.gen) if Int(row[8]) == 1)
  for row in eachrow(mpc.bus)
    bus_i = Int(row[1])   # BUS_I
    btype = Int(row[2])   # BUS_TYPE
    if bus_i in gen_buses && btype != 3
      row[2] = 1.0
    end
  end
  return mpc
end

function run_case(case_name::String, mpc; verbose::Int = 0, show_generator_setup::Bool = false)
  net = with_logger(NullLogger()) do
    Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)
  end
  ite, erg, etime = run_net_acpflow(net = net, max_ite = 30, tol = 1e-8, verbose = verbose, opt_sparse = true, method = :rectangular, show_results = false)
  println("=== $case_name ===")
  println("Converged: ", erg == 0, " (iterations: ", ite, ")")
  if show_generator_setup
    generator_summary(net)
  end
  printACPFlowResults(net, etime, ite, 1e-8; converged = (erg == 0), solver = :rectangular)
  println("")
  return net
end

function run_example_matpower_pv_vs_controller_pq(; verbose::Int = 0, show_generator_setup::Bool = false)
  casefile = "case9.m"
  filename = ensure_case_available!(casefile)
  if isnothing(filename)
    println("Using inline MATPOWER case9 data")
    mpc_pv = case9_inline()
  else
    println("Using MATPOWER case: $filename")
    mpc_pv = Sparlectra.MatpowerIO.read_case(filename; legacy_compat = true)
  end

  net_pv = run_case("Case A: MATPOWER generator buses as PV", mpc_pv; verbose = verbose, show_generator_setup = show_generator_setup)

  mpc_pq = deepcopy(mpc_pv)
  force_generator_buses_to_pq!(mpc_pq)
  net_pq = run_case("Case B: non-slack generators as PQ + constant P(U)/Q(U)", mpc_pq; verbose = verbose, show_generator_setup = show_generator_setup)

  println("=== Delta table (Case B - Case A) at generator buses ===")
  gen_buses = sort(unique(getPosumerBusIndex(ps) for ps in net_pv.prosumpsVec if isGenerator(ps)))
  println(" bus | Vm PV-case | Vm PQ-case | ΔVm     | Va PV-case | Va PQ-case | ΔVa [deg]")
  println("-----+------------+------------+---------+------------+------------+----------")
  for bus in gen_buses
    vm_pv = net_pv.nodeVec[bus]._vm_pu
    vm_pq = net_pq.nodeVec[bus]._vm_pu
    va_pv = net_pv.nodeVec[bus]._va_deg
    va_pq = net_pq.nodeVec[bus]._va_deg
    @printf(" %3d | %10.5f | %10.5f | %+0.5f | %10.4f | %10.4f | %+0.4f\n", bus, vm_pv, vm_pq, vm_pq - vm_pv, va_pv, va_pq, va_pq - va_pv)
  end

  return nothing
end

_ = Base.invokelatest(run_example_matpower_pv_vs_controller_pq)
