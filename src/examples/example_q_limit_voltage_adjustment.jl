using Sparlectra

function build_case(; qmin = -15.0, qmax = 15.0)
  net = Net(name = "q_limit_adjust_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "SLACK", vn_kV = 20.0)
  addBus!(net = net, busName = "PV1", vn_kV = 20.0)
  addBus!(net = net, busName = "PQ1", vn_kV = 20.0)

  addACLine!(net = net, fromBus = "SLACK", toBus = "PV1", length = 5.0, r = 0.1, x = 0.2, c_nf_per_km = 8.5, tanδ = 0.0)
  addACLine!(net = net, fromBus = "PV1", toBus = "PQ1", length = 8.0, r = 0.1, x = 0.2, c_nf_per_km = 8.5, tanδ = 0.0)

  addProsumer!(net = net, busName = "SLACK", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SLACK")
  addProsumer!(net = net, busName = "PV1", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = 1.025, qMin = qmin, qMax = qmax)
  addProsumer!(net = net, busName = "PQ1", type = "ENERGYCONSUMER", p = 70.0, q = 20.0)
  return net
end

function set_controller!(net::Net, bus_name::String; vstep_pu::Union{Nothing,Float64}, tap_steps_down::Union{Nothing,Int}, tap_steps_up::Union{Nothing,Int})
  bus = geNetBusIdx(net = net, busName = bus_name)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    getPosumerBusIndex(ps) == bus || continue
    ps.vstep_pu = vstep_pu
    ps.tap_steps_down = tap_steps_down
    ps.tap_steps_up = tap_steps_up
    return
  end
end

function report_case(label::String, net::Net, mode::Symbol)
  bus = geNetBusIdx(net = net, busName = "PV1")
  node = net.nodeVec[bus]
  status = getNodeType(node) == Sparlectra.PV ? "PV" : "PQ"
  println("$label")
  println("  mode: $mode")
  println("  bus: PV1")
  println("  final voltage pu: $(round(node._vm_pu, digits = 6))")
  println("  final Q MVar: $(round(node._qƩGen, digits = 6))")
  println("  final status: $status")
end

function main()
  # Case A — successful voltage adjustment
  net_a = build_case(qmin = -30.0, qmax = 30.0)
  set_controller!(net_a, "PV1"; vstep_pu = 0.005, tap_steps_down = 20, tap_steps_up = 2)
  runpf!(net_a, 40, 1e-8, 1; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 20)
  report_case("Case A (voltage adjustment succeeds)", net_a, :adjust_vset)

  # Case B — no configuration (immediate PV→PQ fallback)
  net_b = build_case(qmin = -15.0, qmax = 15.0)
  runpf!(net_b, 40, 1e-8, 1; method = :rectangular, qlimit_mode = :adjust_vset)
  report_case("Case B (no controller configured)", net_b, :adjust_vset)

  # Case C — adjustment attempted, no steps available, then PV→PQ fallback
  net_c = build_case(qmin = -15.0, qmax = 15.0)
  set_controller!(net_c, "PV1"; vstep_pu = 0.005, tap_steps_down = 0, tap_steps_up = 0)
  runpf!(net_c, 40, 1e-8, 1; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 1)
  report_case("Case C (controller exhausted -> fallback)", net_c, :adjust_vset)
end

Base.invokelatest(main)
