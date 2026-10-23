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

# file: test/test_legacy_cgmes_importer.jl
# purpose: legacy profile: the testsets of test/test_cgmes_importer.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_cgmes_importer.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_cgmes_importer_tests()
  # cut from testset: CGMES importer
  (function ()
      store = _cgmes_synth_store()
      # The service legs on in-memory deliveries: the synthetic 3W delivery
      # packed as a ZIP is the healthy case, the synthetic EQ with an absent
      # prerequisite the broken one, and a MATPOWER case covers the non-CGMES
      # rejections (the healthy runs on a checked-in delivery are above).
      # --- self-built deliveries (data/cgmes_demo) --------------------------
      #
      # The three checked-in deliveries are exports of shipped cases, written
      # by tools/gen_cgmes_fixtures.jl: sp_case14 (OLTC voltage controller),
      # sp_case118 (54 machines with asymmetric Q limits, synchronous
      # condensers), sp_casePST (phase-shifting transformer, bus link). Their
      # SV profile is the source's solution with every regulated unit holding
      # its target (Q-limits off in the fixture solve), so the profile is the
      # solution of the delivery itself and an import that rebuilds the model
      # reproduces it to floating-point noise. Every solve below keeps
      # Q-limits off for the same reason; what the Q-limit switching makes of
      # the delivery is a solver question, not an importer one.
      solve_fixture!(net; kwargs...) = runpf!(net, 60, 1e-8, 0; method = :rectangular, qlimits_enabled = false, kwargs...)

      @testset "multi-valued references overflow into refsAll (#294 point 9)" begin (function ()
        island = only(objectsOf(store, :TopologicalIsland))
        # refs carries the FIRST value under the short key (plus the dotted-name
        # fallback), which ref() reads, never a silent last-wins
        @test island.refs[:TopologicalNodes] == "_tn_a"
        @test Sparlectra.CGMESImporter.refsAll(store, island, :TopologicalNodes) == ["_tn_a", "_tn_b", "_tn_c"]
        # single-valued references stay single, absent ones come back empty
        line = only(objectsOf(store, :ACLineSegment))
        @test Sparlectra.CGMESImporter.refsAll(store, line, :BaseVoltage) == ["_bv110"]
        @test Sparlectra.CGMESImporter.refsAll(store, line, :NoSuchRef) == String[]
      end)() end

      @testset "ReactiveCapabilityCurve interpolation (#294 point 1)" begin (function ()
        # MicroGrid BE-G1 shape: P −100 → ±200, P 0 → ±300, P 100 → ±200 MVAr.
        pts = [(-100.0, -200.0, 200.0), (0.0, -300.0, 300.0), (100.0, -200.0, 200.0)]
        # interior interpolation (the machine's actual operating point)
        @test Sparlectra.CGMESImporter._curveQHull(pts, -90.0) == (-210.0, 210.0)
        @test Sparlectra.CGMESImporter._curveQHull(pts, 0.0) == (-300.0, 300.0)
        @test Sparlectra.CGMESImporter._curveQHull(pts, 50.0) == (-250.0, 250.0)
        # clamping outside the curve's P domain
        @test Sparlectra.CGMESImporter._curveQHull(pts, -500.0) == (-200.0, 200.0)
        @test Sparlectra.CGMESImporter._curveQHull(pts, 500.0) == (-200.0, 200.0)
        # the sign-convention hull: a curve written in the load convention
        # (y1/y2 swapped in sign) yields the same band
        swapped = [(-100.0, 200.0, -200.0), (0.0, 300.0, -300.0), (100.0, 200.0, -200.0)]
        @test Sparlectra.CGMESImporter._curveQHull(swapped, -90.0) == (-210.0, 210.0)
        # degenerate inputs are refused, the caller falls back to scalars
        @test Sparlectra.CGMESImporter._curveQHull(NTuple{3,Float64}[], 0.0) === nothing
        @test Sparlectra.CGMESImporter._curveQHull([(0.0, 0.0, 0.0)], 0.0) === nothing
        # single-point curve with a real range works
        @test Sparlectra.CGMESImporter._curveQHull([(0.0, -50.0, 50.0)], 25.0) == (-50.0, 50.0)
      end)() end

      # A single ReliCapGrid model is one area of a multi-area system, so its
      # border nodes hang free when imported alone. The combined aliases fetch
      # several areas plus their shared boundary files as ONE delivery. Table
      # consistency is checked here without touching the network or the cache;
      # the fetch itself needs GitHub and is not exercised by the tests.
      @testset "ReliCapGrid combined aliases" begin (function ()
        combined = Sparlectra.CGMESImporter.RELICAPGRID_COMBINED
        singles = Sparlectra.CGMESImporter.RELICAPGRID_ALIASES
        @test haskey(combined, "relicapgrid_cgm")
        @test haskey(combined, "svedala_neighbours")
        # Every member must be a known single model, else the fetch would fail
        # only at download time.
        for (name, members) in combined
          @test !isempty(members)
          for m in members
            @test haskey(singles, m)
          end
        end
        # relicapgrid_cgm spans every single model we offer.
        @test sort(combined["relicapgrid_cgm"]) == sort(collect(keys(singles)))

        # A border file is named Boundary_Border-<AreaA>-<AreaB>.xml. The border is
        # closed by a combination when BOTH its areas are members; otherwise the
        # nodes on that border still hang free.
        border_areas(file) = split(replace(replace(file, "Boundary_Border-" => ""), ".xml" => ""), "-")
        function open_borders(members)
          model_names = Set(lowercase(singles[m].model) for m in members)
          borders = Set{String}()
          for m in members
            union!(borders, singles[m].boundary)
          end
          return sort([b for b in borders if !all(a -> lowercase(a) in model_names, border_areas(b))])
        end

        # The whole point of the CGM alias: not a single border left open.
        @test isempty(open_borders(combined["relicapgrid_cgm"]))
        # The cheap variant closes Svedala's own two borders but inherits its
        # neighbours' outer borders — documented, not accidental.
        @test open_borders(combined["svedala_neighbours"]) == ["Boundary_Border-Espheim-Portheim.xml", "Boundary_Border-Galia-Belgovia.xml"]
        # Svedala alone: both of its borders hang free. That is the situation the
        # combined aliases exist to fix.
        @test length(open_borders(["svedala"])) == 2
      end)() end

      @testset "3W star delivery: leg tap controller (#294 point 4) and star start" begin (function ()
        dir = _cgmes_synth3w_dir()
        # without tap_control: no controllers, taps stay at their SSH position
        plain = importCGMES(path = dir, name = "synth3w_plain")
        @test isempty(collect(Sparlectra._tap_controllers(plain.net)))
        # with tap_control: the star-equivalent leg gets the controller
        res = importCGMES(path = dir, name = "synth3w", tap_control = true)
        @test !any(m -> occursin("not wired in Stage 2", m), res.messages)
        @test any(m -> occursin("tap control: T3W (leg 2)", m), res.messages)
        ctrls = collect(Sparlectra._tap_controllers(res.net))
        @test length(ctrls) == 1
        oltc = only(ctrls)
        @test oltc.mode == :voltage && oltc.control_ratio && oltc.enabled
        # each star leg is addressed by its end (the transformer itself has
        # three branches, so no entry), the tap changer by its controller; the
        # tapped leg runs from its side bus to the star, the others from the
        # star to their side bus
        names = Sparlectra._bus_name_by_idx(res.net)
        for (end_mrid, bus, tapped) in (("_t3w_e1", "HV_BUS", false), ("_t3w_e2", "MV_BUS", true), ("_t3w_e3", "LV_BUS", false))
          kind, k = res.element_by_mrid[end_mrid]
          @test kind == :transformer
          br = res.net.branchVec[k]
          @test names[tapped ? br.fromBus : br.toBus] == bus
          @test startswith(names[tapped ? br.toBus : br.fromBus], "AUX3WT_")
        end
        @test !haskey(res.element_by_mrid, "_t3w")
        @test res.element_by_mrid["_rtc2"] == (:tap_controller, 1)
        @test oltc.trafo == string(res.net.branchVec[last(res.element_by_mrid["_t3w_e2"])].branchIdx)
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
        runres = run_sparlectra(net = res.net, config = cfg)
        @test runres.numerical_converged
        @test latest_control_result(res.net).status == :converged
        @test oltc.converged
        mv = res.net.nodeVec[res.net.busDict["MV_BUS"]]
        # target 1.04 pu ± 0.0118 pu — and the tap really moved: the
        # uncontrolled operating point sits at ≈1.0135 pu, outside the band
        @test abs(mv._vm_pu - 1.04) <= 2.6 / 110.0 / 2.0 + 1e-9
        @test mv._vm_pu > 1.025

        # Star start: the star point has no SvVoltage and starts at the value
        # that closes its own node balance against its sides. An end-1 tap 8 %
        # off neutral puts that value far from the HV voltage, the provisional
        # start. The SV is written from the solved model, so at
        # the start the star carries no current and the three leg flows of
        # compareWithSV equal the SV. Every end carries a magnetizing
        # admittance, so the star balance includes the arms at the star.
        tdir = _cgmes_synth3w_dir(end1_tap_step = 19, magnetizing = true)
        model = importCGMES(path = tdir, name = "synth3w_tap1_model")
        # the unregulated end-1 changer (step 19 of 1..21, 1 % per step,
        # neutral 11) gives its leg the linear range around the nominal ratio,
        # declared as the linear grid the tap estimation fixes on
        leg1 = model.net.branchVec[last(model.element_by_mrid["_t3w_e1"])]
        rho = leg1.tap_ratio / 1.08
        @test all(isapprox.((leg1.tap_min, leg1.tap_max, leg1.tap_step), (0.9 * rho, 1.1 * rho, 0.01 * rho); rtol = 1e-12)) && leg1.tap_ratio_points === nothing
        @test leg1.tap_ratio_grid === :linear
        # the magnetizing admittance of each end sits at the star node (0.42.0),
        # for the tapped legs 1 and 2 (side bus -> star) and the untapped leg 3
        # (star -> side bus): the arm seen at the star is the end's g + jb
        # referred to the star base, the arm at the side bus is zero
        mstar = model.net.busDict["AUX3WT_T3W"]
        for (k, (Uk, (g, b))) in enumerate(zip((220.0, 110.0, 20.0), _CGMES_SYNTH3W_MAG))
          br = model.net.branchVec[last(model.element_by_mrid["_t3w_e$(k)"])]
          y0f = Sparlectra._branch_y0_from(br) / abs2(Sparlectra.calcBranchRatio(br))
          y0t = Sparlectra._branch_y0_to(br)
          at_star, at_side = br.fromBus == mstar ? (y0f, y0t) : (y0t, y0f)
          @test (br.fromBus == mstar) == (k == 3)
          @test at_star ≈ complex(g, b) * Uk^2 / model.net.baseMVA rtol = 1e-12
          @test iszero(at_side)
        end
        @test runpf!(model.net, 30, 1e-10, 0)[2] == 0
        Vm = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in model.net.nodeVec]
        sv_v = Dict(tn => (abs(Vm[model.net.busDict[bus]]) * vn, rad2deg(angle(Vm[model.net.busDict[bus]]))) for (bus, tn, vn) in (("HV_BUS", "_tn_hv", 220.0), ("MV_BUS", "_tn_mv", 110.0), ("LV_BUS", "_tn_lv", 20.0)))
        sv_f = Dict{String,Tuple{Float64,Float64}}()
        for k in 1:3
          br = model.net.branchVec[last(model.element_by_mrid["_t3w_e$(k)"])]
          # the flow at the side bus: a tapped leg (ends 1 and 2 here) runs side
          # bus -> star, so its terminal is the from side
          side = last(model.branch_side_of_terminal["_t_t3w_$(k)"])
          s = Sparlectra._sv_branch_terminal_flow(br, Vm, side, model.net.baseMVA)
          sv_f["_t_t3w_$(k)"] = (real(s), imag(s))
        end
        write(joinpath(tdir, "synth3w_SV.xml"), _cgmes_synth3w_sv(sv_v, sv_f))
        res1 = importCGMES(path = tdir, name = "synth3w_tap1")
        star = res1.net.busDict["AUX3WT_T3W"]
        V0 = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in res1.net.nodeVec]
        @test abs((Sparlectra.createYBUS(net = res1.net, sparse = true) * V0)[star]) < 1e-9
        legs = filter(r -> r.kind in (:branch_from, :branch_to), compareWithSV(res1).flows.rows)
        @test length(legs) == 3
        @test all(r -> abs(r.dp) < 1e-6 && abs(r.dq) < 1e-6, legs)
        # a side without a usable SV voltage (v = 0) keeps the end-1 SV start
        sv_v["_tn_lv"] = (0.0, 0.0)
        write(joinpath(tdir, "synth3w_SV.xml"), _cgmes_synth3w_sv(sv_v, sv_f))
        res0 = importCGMES(path = tdir, name = "synth3w_tap1_nolv")
        nd0 = res0.net.nodeVec[res0.net.busDict["AUX3WT_T3W"]]
        @test nd0._vm_pu == sv_v["_tn_hv"][1] / 220.0 && nd0._va_deg == sv_v["_tn_hv"][2]

        # An open side terminal (Terminal.connected = false, the unit in
        # service) keeps its leg open at the side bus: the winding is still
        # magnetized from the star, so the end's admittance stays at the star
        # as the input admittance of the open leg. End 2 is tapped (side bus ->
        # star, the side is the from end), end 3 untapped (star -> side bus).
        # The side bus has no other feed and is de-energized.
        for k in (2, 3)
          odir = _cgmes_synth3w_dir(end1_tap_step = 19, magnetizing = true, open_end = k)
          ores = importCGMES(path = odir, name = "synth3w_open$(k)")
          onet = ores.net
          ostar = onet.busDict["AUX3WT_T3W"]
          obr = onet.branchVec[last(ores.element_by_mrid["_t3w_e$(k)"])]
          @test obr.status == 0
          @test Sparlectra._branch_terminal_state(obr) == (k == 3 ? :open_to : :open_from)
          @test (obr.fromBus == ostar) == (k == 3)
          g, b = _CGMES_SYNTH3W_MAG[k]
          ymag = complex(g, b) * (220.0, 110.0, 20.0)[k]^2 / onet.baseMVA
          @test Sparlectra._open_terminal_yin(obr) ≈ ymag rtol = 1e-9
          @test runpf!(onet, 30, 1e-10, 0)[2] == 0
          # star balance: the two closed legs feed the open leg's magnetizing
          # admittance at the star voltage
          V = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in onet.nodeVec]
          i_legs = zero(ComplexF64)
          for j in 1:3
            j == k && continue
            br = onet.branchVec[last(ores.element_by_mrid["_t3w_e$(j)"])]
            y11, y12, y21, y22 = Sparlectra.calcAdmittance(br, 1.0, 1.0)
            i_legs += br.fromBus == ostar ? y11 * V[br.fromBus] + y12 * V[br.toBus] : y21 * V[br.fromBus] + y22 * V[br.toBus]
          end
          @test abs(i_legs + ymag * V[ostar]) < 1e-9
          # export round trip: the open terminal goes out as connected = false
          # and comes back as the leg open at its side bus. The export flattens
          # the taps into the end ratedU and scales the star voltage by the
          # end-1 ratio, so the re-imported leg runs star -> side bus and its
          # arm differs in per unit; the power it draws at the star does not.
          files = writeCGMESFiles(onet; path = mktempdir(), notices = String[])
          back = importCGMES(path = dirname(files[1]), name = "synth3w_open$(k)_back")
          @test !Sparlectra.CGMESImporter.terminalConnected(back.store.objects["_t_t3w_$(k)"])
          bbr = back.net.branchVec[last(back.element_by_mrid["_t3w_e$(k)"])]
          bst = Sparlectra._branch_terminal_state(bbr)
          @test bst in (:open_from, :open_to)
          @test (bst == :open_to ? bbr.toBus : bbr.fromBus) == back.net.busDict[k == 2 ? "MV_BUS" : "LV_BUS"]
          @test runpf!(back.net, 30, 1e-10, 0)[2] == 0
          Vb = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in back.net.nodeVec]
          bstar = back.net.busDict["AUX3WT_T3W"]
          @test abs2(Vb[bstar]) * Sparlectra._open_terminal_yin(bbr) ≈ abs2(V[ostar]) * ymag rtol = 1e-8
          for bus in ("HV_BUS", k == 2 ? "LV_BUS" : "MV_BUS")
            n1 = onet.nodeVec[onet.busDict[bus]]
            n2 = back.net.nodeVec[back.net.busDict[bus]]
            @test abs(n1._vm_pu - n2._vm_pu) < 1e-9 && abs(n1._va_deg - n2._va_deg) < 1e-7
          end
        end
      end)() end

      @testset "Remote-regulating machine (#294 point 3)" begin (function ()
        # default: the remote control stays the Stage-1 held-PV fallback
        plain = importCGMES(path = _cgmes_synth_rvc_dir(), name = "synthrvc_plain")
        @test any(m -> occursin("G_B has a remote voltage RegulatingControl", m) && occursin("held PV at its own bus", m), plain.messages)
        @test isempty(Sparlectra._machine_controllers(plain.net))
        b_plain = plain.net.nodeVec[plain.net.busDict["BUS_B"]]
        @test b_plain._nodeType == Sparlectra.PV

        # machine_control = true: G_B becomes a PQ injection with an outer-loop
        # controller regulating BUS_C to 112.2 kV = 1.02 pu
        res = importCGMES(path = _cgmes_synth_rvc_dir(), name = "synthrvc", machine_control = true)
        @test any(m -> occursin("machine control: G_B", m) && occursin("1.02 pu at BUS_C", m), res.messages)
        ctrls = Sparlectra._machine_controllers(res.net)
        @test length(ctrls) == 1
        ctrl = only(ctrls)
        @test ctrl.target_bus == "BUS_C"
        @test isapprox(ctrl.target_vm_pu, 112.2 / 110.0; atol = 1e-9)
        @test ctrl.qmin_mvar == -50.0 && ctrl.qmax_mvar == 50.0
        @test res.net.nodeVec[res.net.busDict["BUS_B"]]._nodeType == Sparlectra.PQ
        # the control loop drives the remote bus onto the target
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
        runres = run_sparlectra(net = res.net, config = cfg)
        @test runres.numerical_converged
        @test latest_control_result(res.net).status == :converged
        @test ctrl.converged
        @test abs(res.net.nodeVec[res.net.busDict["BUS_C"]]._vm_pu - 112.2 / 110.0) <= ctrl.deadband_vm_pu + 1e-9

        # fallback: with G_C reconnected the target bus is already voltage-held
        # locally — the plan reverts to held-PV with a notice
        held = importCGMES(path = _cgmes_synth_rvc_dir(sm2_connected = true), name = "synthrvc_held", machine_control = true)
        @test any(m -> occursin("G_B", m) && occursin("remote voltage control not attachable", m) && occursin("already voltage-held", m), held.messages)
        @test isempty(Sparlectra._machine_controllers(held.net))
        @test held.net.nodeVec[held.net.busDict["BUS_B"]]._nodeType == Sparlectra.PV

        # controlEnabled = false on G_B: no controller, no held-PV fallback, a
        # PQ unit at its SSH operating point
        off_store = loadCGMES(_cgmes_synth_rvc_dir())
        off_store.objects["_smb"].attrs[:controlEnabled] = "false"
        off = importCGMES(path = off_store, name = "synthrvc_off", machine_control = true)
        @test isempty(Sparlectra._machine_controllers(off.net))
        @test off.net.nodeVec[off.net.busDict["BUS_B"]]._nodeType == Sparlectra.PQ
        @test !any(m -> occursin("G_B", m) && (occursin("machine control", m) || occursin("remote voltage", m)), off.messages)
        g_b = off.net.prosumpsVec[last(off.element_by_mrid["_smb"])]
        @test (g_b.pVal, g_b.qVal, g_b.minQ, g_b.maxQ) == (20.0, 0.0, -50.0, 50.0)
        # with machine_control = true G_B was a PQ injection before the flag was
        # read; the default import is where the flag changes the bus type: no
        # held-PV fallback, BUS_B is PQ
        off_plain = importCGMES(path = off_store, name = "synthrvc_off_plain")
        @test off_plain.net.nodeVec[off_plain.net.busDict["BUS_B"]]._nodeType == Sparlectra.PQ
        @test !Sparlectra.isRegulating(off_plain.net.prosumpsVec[last(off_plain.element_by_mrid["_smb"])])
        @test !any(m -> occursin("G_B has a remote voltage RegulatingControl", m), off_plain.messages)
      end)() end

      @testset "NonlinearShuntCompensator mapping (#294 point 7)" begin (function ()
        dir = _cgmes_synth3w_dir()
        res = importCGMES(path = dir, name = "synth3w_nlsh")
        # sections = 2 of 3 points: B = (1.0 + 0.8)e-4 S * 110² = 2.178 MVAr
        @test any(m -> occursin("nonlinear shunt NLSH", m) && occursin("sections 2/3", m) && occursin("2.178 MVAr", m), res.messages)
        @test !any(m -> occursin("skip: NonlinearShuntCompensator", m), res.messages)
        # multi-valued reference guard (#294 point 9): the SV TopologicalIsland
        # membership list triggers exactly one class/property notice
        @test count(m -> occursin("multi-valued reference TopologicalIsland.TopologicalNodes on 1 object(s)", m), res.messages) == 1
        @test length(res.net.shuntVec) == 1
        @test res.element_by_mrid["_nlsh"] == (:shunt, 1)
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        # the SV comparison has a flow row for it: an SvPowerFlow added to a
        # copy of the store at the model value for the SV voltage of MV_BUS
        # (load sign, capacitive Q negative) is met at the SV start
        local nl_store = loadCGMES(dir)
        local q_sv = -2.178 * (111.5 / 110.0)^2
        nl_store.objects["_nlsh_svpf"] = Sparlectra.CGMESImporter.CIMObject("_nlsh_svpf", :SvPowerFlow, Dict{Symbol,String}(:p => "0", :q => string(q_sv)), Dict{Symbol,String}(:Terminal => "_t_nlsh"), "nlsh_test")
        push!(get!(Vector{String}, nl_store.byclass, :SvPowerFlow), "_nlsh_svpf")
        local nl_row = only(r for r in compareWithSV(importCGMES(path = nl_store, name = "synth3w_nlsh_sv")).flows.rows if r.kind == :shunt)
        @test nl_row.name == "NLSH"
        @test abs(nl_row.dp) < 1e-9 && abs(nl_row.dq) < 1e-9
      end)() end

      # Stage 2 semantics on the exporter's TapChangerControl: without
      # tap_control the delivery imports as a fixed-tap model, with it the
      # OLTC of sp_case14 comes back as the voltage controller the source case
      # carries (same target bus, target and deadband), and the control loop
      # meets the target within its deadband.
      @testset "OLTC controller of sp_case14 arrives as a tap controller (Stage 2)" begin (function ()
        dir = cgmes_fixture_dir("sp_case14")
        stage1 = importCGMES(path = dir, name = "s1")
        @test isempty(collect(Sparlectra._tap_controllers(stage1.net)))
        res = importCGMES(path = dir, tap_control = true, name = "s2")
        ctrls = collect(Sparlectra._tap_controllers(res.net))
        @test length(ctrls) == 1
        c = only(ctrls)
        @test c.mode == :voltage
        @test c.control_ratio
        @test !c.control_phase
        @test c.enabled
        @test isempty(c.followers)
        @test any(m -> startswith(m, "tap control:") && occursin("voltage", m), res.messages)
        src = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
        sc = only(collect(Sparlectra._tap_controllers(src)))
        @test c.target_bus == sc.target_bus
        @test isapprox(something(c.target_vm_pu, NaN), sc.target_vm_pu; atol = 1e-9)
        @test isapprox(c.deadband_vm_pu, sc.deadband_vm_pu; atol = 1e-9)
        # the controlled transformer connects the same buses as in the source
        # and keeps its tap machinery (range and step from the RatioTapChanger)
        names = Sparlectra._bus_name_by_idx(res.net)
        src_names = Sparlectra._bus_name_by_idx(src)
        br = Sparlectra._find_trafo_branch(res.net, c.trafo)
        sbr = Sparlectra._find_trafo_branch(src, sc.trafo)
        @test Set([names[br.fromBus], names[br.toBus]]) == Set([src_names[sbr.fromBus], src_names[sbr.toBus]])
        @test br.has_ratio_tap
        @test br.tap_step > 0.0
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
        runres = run_sparlectra(net = res.net, config = cfg)
        @test runres.numerical_converged
        @test latest_control_result(res.net).status == :converged
        @test c.converged
        @test abs(c.achieved_vm_pu - c.target_vm_pu) <= c.deadband_vm_pu
      end)() end

      # The PST of sp_casePST travels as a single-step PhaseTapChangerLinear on
      # end 1; its shift comes back with the source sign, and the closed
      # Breaker comes back as the bus link.
      @testset "PST of sp_casePST keeps its tap-side sign, bus link survives" begin (function ()
        res = importCGMES(path = cgmes_fixture_dir("sp_casePST"), name = "pst")
        src = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
        shifted(net) = [br for br in net.branchVec if br.phase_shift_deg != 0.0]
        @test length(shifted(src)) == 1
        @test length(shifted(res.net)) == 1
        sbr, ibr = only(shifted(src)), only(shifted(res.net))
        @test sbr.phase_shift_deg < 0.0
        @test isapprox(ibr.phase_shift_deg, sbr.phase_shift_deg; atol = 1e-9)
        @test isapprox(ibr.angle, sbr.angle; atol = 1e-9)
        names = Sparlectra._bus_name_by_idx(res.net)
        src_names = Sparlectra._bus_name_by_idx(src)
        @test (names[ibr.fromBus], names[ibr.toBus]) == (src_names[sbr.fromBus], src_names[sbr.toBus])
        @test any(m -> occursin("bus link", m) && occursin("Breaker", m), res.messages)
        @test length(res.net.linkVec) == 1
        @test only(res.net.linkVec).status == 1
        # the link is addressed under the Breaker's mRID
        links = [mrid for (mrid, (kind, _)) in res.element_by_mrid if kind == :link]
        @test length(links) == 1
        @test res.store.objects[only(links)].class == :Breaker
        @test res.element_by_mrid[only(links)] == (:link, 1)
      end)() end

      # SCF export of a delivery (issue #342): a delivery-sourced net reaches
      # the case format with its source identity intact, the file is
      # deterministic, and the case re-imports to the same solution.
      @testset "SCF export of a fixture delivery round-trips" begin (function ()
        res = importCGMES(path = cgmes_fixture_dir("sp_case14"), name = "scf14")
        @test solve_fixture!(res.net)[2] == 0
        root = Sparlectra.net_to_scf(res.net; source_format = "cgmes", source_reference = "sp_case14 fixture")
        extra = root["sparlectra"]["extra"]
        tagged = [e for e in values(extra) if get(e, "source_id_kind", "") == "cgmes_mrid"]
        @test !isempty(tagged)
        @test all(e -> !startswith(String(e["external_id"]), "#"), tagged)
        names = Set(String(e["name"]) for e in values(extra))
        @test all(bus -> bus in names, keys(res.net.busDict))
        d = mktempdir()
        a = exportSCF(res.net; file = joinpath(d, "a.scf.json"), source_format = "cgmes")
        b = exportSCF(res.net; file = joinpath(d, "b.scf.json"), source_format = "cgmes")
        @test read(a, String) == read(b, String)
        back = importSCF(a)
        @test length(back.nodeVec) == 14
        @test solve_fixture!(back)[2] == 0
        for (bus, i) in res.net.busDict
          j = back.busDict[bus]
          @test isapprox(res.net.nodeVec[i]._vm_pu, back.nodeVec[j]._vm_pu; atol = 1e-6)
          @test isapprox(res.net.nodeVec[i]._va_deg, back.nodeVec[j]._va_deg; atol = 1e-6)
        end

        # A CGMES changer steps linearly from its import ratio (tap_ratio_grid
        # :linear). The file carries the position on that grid (`grid:
        # "linear"`), so a tap a controller moved off the import ratio comes
        # back exactly. Written on the fraction grid, 8 steps of 1 % from the
        # neutral import ratio came back as 1.0753 instead of 1.08 (rounded
        # -7 fraction steps), the band 0.9/1.1 as 0.9009/1.0989, and the
        # voltages moved by 4.2e-3 pu.
        t3 = importCGMES(path = _cgmes_synth3w_dir(end1_tap_step = 11), name = "scf3w_lin")
        leg = t3.net.branchVec[last(t3.element_by_mrid["_t3w_e1"])]
        @test leg.tap_ratio_grid === :linear
        leg.tap_ratio = Sparlectra._ratio_grid_shift(leg, leg.ratio, 8)
        @test leg.tap_ratio == leg.ratio + 8 * leg.tap_step
        @test solve_fixture!(t3.net)[2] == 0
        f3 = exportSCF(t3.net; file = joinpath(d, "lin.scf.json"), source_format = "cgmes")
        root3 = Sparlectra.net_to_scf(t3.net; source_format = "cgmes")
        leg_id = parse(Int, only(k for (k, e) in root3["sparlectra"]["extra"] if get(e, "branch_index", 0) == leg.branchIdx))
        row3 = only(r for r in root3["sparlectra"]["components"]["tap_changer"] if r["branch"] == leg_id)
        c1 = only(c for c in row3["controllers"] if c["index"] == 1)
        @test (c1["pos"], c1["pos_min"], c1["pos_max"], c1["grid"]) == (8, -10, 10, "linear")
        back3 = importSCF(f3)
        bleg = back3.branchVec[findfirst(b -> getCompName(b.comp) == getCompName(leg.comp), back3.branchVec)]
        @test bleg.tap_ratio_grid === :linear
        @test (bleg.tap_ratio, bleg.tap_min, bleg.tap_max, bleg.tap_step) == (leg.tap_ratio, leg.tap_min, leg.tap_max, leg.tap_step)
        @test solve_fixture!(back3)[2] == 0
        for (bus, i) in t3.net.busDict
          j = back3.busDict[bus]
          @test isapprox(t3.net.nodeVec[i]._vm_pu, back3.nodeVec[j]._vm_pu; atol = 1e-12)
        end
        # the re-read case is a fixed point of the writer
        @test read(exportSCF(back3; file = joinpath(d, "lin2.scf.json"), source_format = "cgmes"), String) == read(exportSCF(importSCF(joinpath(d, "lin2.scf.json")); file = joinpath(d, "lin3.scf.json"), source_format = "cgmes"), String)
      end)() end

      # The framework/API path: a delivery zip as casefile, dispatch by
      # auto-detection, cgmes_import_report.md next to run.log, and the SV comparison
      # artifacts of the default start (auto resolves to sv on a delivery
      # that carries SvVoltage).
      @testset "run_sparlectra_api dispatch on a fixture delivery zip" begin (function ()
        z = _pack_cgmes_fixture_zip("sp_case14")
        out = mktempdir()
        r = run_sparlectra_api(casefile = z, output_dir = out)
        @test r.status == :succeeded
        @test r.metadata["input_format_detected"] == "cgmes"
        @test r.metadata["cgmes_version"] == "2.4.15"
        @test r.metadata["cgmes_buses"] == 14
        @test r.metadata["cgmes_slack_bus"] == "Ostheim_110"
        @test r.metadata["cgmes_messages"] > 0
        @test isfile(joinpath(out, "cgmes_import_report.md"))
        @test !isfile(joinpath(out, "cgmes.log"))
        # the report is a Markdown artifact: the Web UI renders it
        art = only(a for a in collect_sparlectra_api_artifacts(out) if a.name == "cgmes_import_report.md")
        @test art.mime_type == "text/markdown"
        @test r.metadata["cgmes_start_values"] == "sv"
        @test r.metadata["cgmes_sv_compare_status"] == "converged"
        @test r.metadata["cgmes_sv_compare_n"] == 14
        @test r.metadata["cgmes_sv_compare_max_dvm"] < 1e-6
        @test isfile(joinpath(out, "sv_compare.csv"))
        @test isfile(joinpath(out, "sv_compare_flows.csv"))
        @test occursin("CGMES start values: sv (auto: delivery carries SvVoltage for", read(joinpath(out, "run.log"), String))

        # cgmes_import.start_values selects the start state and wins over a
        # hostile power_flow.flatstart; both starts land on the SV solution
        for (mode, hostile_flatstart) in (("flat", "false"), ("sv", "true"))
          out_m = mktempdir()
          cfg_m = joinpath(out_m, "c.yaml")
          write(cfg_m, "config_version: 1\npower_flow:\n  flatstart: " * hostile_flatstart * "\n  tol: 1.0e-10\ncgmes_import:\n  start_values: " * mode * "\n")
          rm_ = run_sparlectra_api(casefile = z, config_file = cfg_m, output_dir = out_m, case_format = :cgmes)
          @test rm_.status == :succeeded
          run_log = read(joinpath(out_m, "run.log"), String)
          @test occursin("CGMES start values: " * mode, run_log)
          @test occursin("overrides: power_flow.flatstart=" * hostile_flatstart, run_log)
          @test rm_.metadata["cgmes_start_values"] == mode
          @test rm_.metadata["cgmes_sv_compare_max_dvm"] < 1e-6
          sv_rows = run_csv(joinpath(out_m, "sv_compare.csv"))
          @test sv_rows.header == ["bus", "vm_pu", "sv_vm_pu", "dvm", "va_deg", "sv_va_deg", "dva", "dva_aligned"]
          @test length(sv_rows.rows) == 14
        end

        # under auto, power_flow.flatstart asks for the flat start although the
        # delivery carries SvVoltage (the Settings page offers that one switch)
        out_a = mktempdir()
        cfg_a = joinpath(out_a, "c.yaml")
        write(cfg_a, "config_version: 1\npower_flow:\n  flatstart: true\n")
        ra = run_sparlectra_api(casefile = z, config_file = cfg_a, output_dir = out_a, case_format = :cgmes)
        @test ra.status == :succeeded
        @test ra.metadata["cgmes_start_values"] == "flat"
        @test occursin("CGMES start values: flat (auto: power_flow.flatstart asks for the flat start)", read(joinpath(out_a, "run.log"), String))

        # the report survives a failed solve: it is written right after the
        # import, not after the power flow (one flat-start iteration cannot
        # reach 1e-14; rescue off, because the rescue ladder's alternate start
        # is the delivery's SV, an exact solution that would converge at once
        # since the projection keeps a seed closer than its DC start, 0.20.4)
        out_f = mktempdir()
        cfg_f = joinpath(out_f, "c.yaml")
        write(cfg_f, "config_version: 1\npower_flow:\n  max_iter: 1\n  tol: 1.0e-14\n  rescue: false\ncgmes_import:\n  start_values: flat\n")
        rf = run_sparlectra_api(casefile = z, config_file = cfg_f, output_dir = out_f)
        @test rf.status != :succeeded
        @test startswith(cgmes_report_plain(out_f), "CGMES import report\n")
      end)() end

      # The Web UI "Analyze import" and "Short circuit" buttons on the
      # fixtures, which carry what their source cases carry. sp_case14 comes
      # from a case file with the feeder of its external grid: the delivery
      # has the feeder's short-circuit currents, the two machines have no
      # data and are flagged. sp_casePST comes from a case without any source
      # data: nothing an impedance could come from, the button is not offered
      # and the run refuses. With a rated power on its machines the table is
      # complete on the default reactance and the run ends with a warning.
      @testset "import_analysis_mode and short_circuit_mode service runs on a fixture" begin (function ()
        root = mktempdir()
        cfg = joinpath(root, "c.yaml")
        write(cfg, "config_version: 1\npower_flow:\n  max_iter: 40\n")
        z = _pack_cgmes_fixture_zip("sp_case14")

        ia = start_powerflow_run(Dict("casefile" => z, "config_file" => cfg, "output_root" => root, "import_analysis_mode" => true))
        @test ia["success"] === true
        @test ia["metadata"]["run_mode"] == "import_analysis"
        @test ia["metadata"]["import_analysis_missing_dependencies"] == 0
        report = read(joinpath(root, ia["run_id"], "import_analysis.txt"), String)
        @test occursin("Supplied models:", report)
        @test occursin("Verdict:", report)

        @test !SparlectraApp._webui_case_has_short_circuit_data(cgmes_fixture_dir("sp_casePST"))
        refused = start_powerflow_run(Dict("casefile" => _pack_cgmes_fixture_zip("sp_casePST"), "config_file" => cfg, "output_root" => root, "short_circuit_mode" => true))
        @test refused["status"] == "failed"
        @test refused["reason"] == "short_circuit_data_missing"
        @test occursin("rated power", refused["message"])
        rated_dir = mktempdir()
        for f in readdir(cgmes_fixture_dir("sp_casePST"); join = true)
          text = read(f, String)
          occursin("_EQ", basename(f)) && (text = replace(text, "  </cim:SynchronousMachine>" => "    <cim:RotatingMachine.ratedS>100</cim:RotatingMachine.ratedS>\n  </cim:SynchronousMachine>"))
          write(joinpath(rated_dir, basename(f)), text)
        end
        @test SparlectraApp._webui_case_has_short_circuit_data(rated_dir)
        rated_zip = joinpath(mktempdir(), "sp_casePST_rated_CGMES.zip")
        ZipArchives.ZipWriter(rated_zip) do w
          for f in sort(readdir(rated_dir; join = true))
            ZipArchives.zip_newfile(w, basename(f))
            write(w, read(f))
          end
        end
        # the engine's warnings about substituted defaults belong to the
        # result: the service writes them into run.log, none reaches the
        # console (an empty list of expected warnings fails on any)
        on_defaults = run_with_expected_warnings(() -> start_powerflow_run(Dict("casefile" => rated_zip, "config_file" => cfg, "output_root" => root, "short_circuit_mode" => true)), String[])
        on_defaults_log = read(joinpath(on_defaults["output_dir"], "run.log"), String)
        @test occursin("Substituted defaults:", on_defaults_log)
        @test occursin("has no usable x''_d", on_defaults_log)
        @test on_defaults["success"] === true
        @test on_defaults["status"] == "warning"
        @test on_defaults["reason"] == "short_circuit_defaults_only"
        # sp_case14: the feeder of the source case is in the delivery
        @test SparlectraApp._webui_case_has_short_circuit_data(cgmes_fixture_dir("sp_case14"))
        @test occursin("ExternalNetworkInjection.maxInitialSymShCCurrent", read(joinpath(cgmes_fixture_dir("sp_case14"), "sp_case14_EQ.xml"), String))
        sc = run_with_expected_warnings(() -> start_powerflow_run(Dict("casefile" => z, "config_file" => cfg, "output_root" => root, "short_circuit_mode" => true)), String[])
        sc_log = read(joinpath(sc["output_dir"], "run.log"), String)
        # the delivery carries the maximum active power of its units
        # (GeneratingUnit.maxOperatingP), the default stands on it
        @test all(text -> occursin(text, sc_log), ("Substituted defaults:", "default 0.2 pu on its maximum active power 50.0 MW", "default 0.2 pu on its maximum active power 35.0 MW"))
        @test sc["success"] === true
        # one source with data, two without: a warning that says so
        @test sc["status"] == "warning"
        @test sc["reason"] == "short_circuit_partial_defaults"
        @test sc["metadata"]["sc_sources"] == 3 && sc["metadata"]["sc_sources_with_data"] == 1
        @test sc["metadata"]["run_mode"] == "short_circuit"
        @test sc["metadata"]["sc_case_rows"] == 14
        @test sc["metadata"]["sc_flagged_rows"] == 14
        @test sc["metadata"]["sc_max_ik_kA"] > 0.0
        rid = sc["run_id"]
        max_rows = run_csv(joinpath(root, rid, "short_circuit_max.csv"))
        @test max_rows.header == ["bus", "vn_kV", "island", "status", "c", "zk_ohm", "rx_ratio", "ik_kA", "sk_MVA", "kappa", "ip_kA", "flagged", "reasons", "lower_bound"]
        @test length(max_rows.rows) == 14
        @test all(occursin("true", row[12]) for row in max_rows.rows)
        # the delivery and the case file it was written from give the same
        # currents: the reference values of the case file hold here as well
        case_reference = Sparlectra.scf_json_parse(read(joinpath(@__DIR__, "fixtures", "demo_cases", "sp_case14.json"), String))["short_circuit"]
        for (bus, reference) in case_reference
          row = only(r for r in max_rows.rows if first(r) == bus)
          @test isapprox(max_rows.number(row[8]), reference["ik_kA"]; atol = 1e-4)
        end
        @test isfile(joinpath(root, rid, "short_circuit_min.csv"))
        run_log = read(joinpath(root, rid, "run.log"), String)
        @test occursin("Short-circuit run", run_log)
        @test occursin("SynchronousMachine", run_log)
      end)() end

  end)()
  return true
end
