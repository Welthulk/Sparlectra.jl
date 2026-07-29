# file: test/test_cgmes_importer.jl
# Tests for the CGMES importer (task_plan_cgmes_import.md, phases C–E):
# generic reader semantics on a synthetic in-memory fixture, plus
# summarizeCGMES assertions on the ENTSO-E MicroGrid when the local test-set
# cache exists (fetched by examples/experimental/cgmes_fetch_testsets.jl;
# skipped otherwise — no network access from tests).

using Sparlectra.CGMESImporter:
  CGMESFile, CGMESStore, collectCGMESFiles, loadCGMES, summarizeCGMES, objectsOf, countOf, num, str, boolval, enumval, ref, unresolvedReferences, readCGMESFile!, CIMObject

const _CGMES_SYNTH_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:eq-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
<md:Model.profile>http://entsoe.eu/CIM/EquipmentShortCircuit/3/1</md:Model.profile>
</md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110">
<cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage>
</cim:BaseVoltage>
<cim:ACLineSegment rdf:ID="_line1">
<cim:IdentifiedObject.name>L1</cim:IdentifiedObject.name>
<cim:ACLineSegment.r>2.5</cim:ACLineSegment.r>
<cim:ACLineSegment.x>10.0</cim:ACLineSegment.x>
<cim:ACLineSegment.gch>0.001</cim:ACLineSegment.gch>
<cim:ACLineSegment.r0>7.5</cim:ACLineSegment.r0>
<cim:Conductor.length>12.0</cim:Conductor.length>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/>
</cim:ACLineSegment>
<cim:SynchronousMachine rdf:ID="_sm1">
<cim:IdentifiedObject.name>G1</cim:IdentifiedObject.name>
<cim:SynchronousMachine.operatingMode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#SynchronousMachineOperatingMode.generator"/>
<cim:RotatingMachine.ratedS>150</cim:RotatingMachine.ratedS>
<cim:SynchronousMachine.earthing>false</cim:SynchronousMachine.earthing>
<cim:Equipment.EquipmentContainer rdf:resource="#_missing_container"/>
</cim:SynchronousMachine>
</rdf:RDF>
"""

const _CGMES_SYNTH_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:ssh-model">
<md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile>
</md:FullModel>
<cim:SynchronousMachine rdf:about="#_sm1">
<cim:RotatingMachine.p>-120</cim:RotatingMachine.p>
<cim:RotatingMachine.q>-30</cim:RotatingMachine.q>
</cim:SynchronousMachine>
</rdf:RDF>
"""

const _CGMES_SYNTH_DIFF = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:dm="http://iec.ch/TC57/61970-552/DifferenceModel/1#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<dm:DifferenceModel rdf:about="urn:uuid:diff-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
</dm:DifferenceModel>
<cim:ACLineSegment rdf:ID="_should_not_appear">
<cim:ACLineSegment.r>1.0</cim:ACLineSegment.r>
</cim:ACLineSegment>
</rdf:RDF>
"""

function _cgmes_synth_store()::CGMESStore
  dir = mktempdir()
  write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
  write(joinpath(dir, "synth_SSH.xml"), _CGMES_SYNTH_SSH)
  write(joinpath(dir, "synth_DIFF.xml"), _CGMES_SYNTH_DIFF)
  return loadCGMES(dir)
end

function run_cgmes_importer_tests()
  @testset "CGMES importer" begin
    store = _cgmes_synth_store()

    @testset "profile classification and version" begin
      @test store.version == "2.4.15"
      eqinfo = only(filter(f -> occursin("EQ", f.name), store.files))
      @test eqinfo.header == :FullModel
      @test eqinfo.profiles == Set([:EQ, :EQ_SC])   # per-file profile *set* (A2/D-7)
      @test !eqinfo.skipped
    end

    @testset "DifferenceModel files are skipped with reason (D-6)" begin
      diffinfo = only(filter(f -> occursin("DIFF", f.name), store.files))
      @test diffinfo.skipped
      @test diffinfo.header == :DifferenceModel
      @test !isempty(diffinfo.skip_reason)
      @test !haskey(store.objects, "_should_not_appear")
    end

    @testset "rdf:ID creates, literals and inherited attributes" begin
      @test countOf(store, :ACLineSegment) == 1
      line = only(objectsOf(store, :ACLineSegment))
      @test str(line, :name) == "L1"                 # IdentifiedObject.name → :name
      @test num(line, :r) == 2.5
      @test num(line, :gch) == 0.001                 # CGMES conductance present
      @test num(line, :r0) == 7.5                    # short-circuit attribute read (§7.7)
      @test num(line, :length) == 12.0               # Conductor.length → :length
      @test num(line, :missing_attr, 99.0) == 99.0
    end

    @testset "references and enums" begin
      line = only(objectsOf(store, :ACLineSegment))
      bv = ref(store, line, :BaseVoltage)
      @test bv !== nothing && bv.class == :BaseVoltage
      @test num(bv, :nominalVoltage) == 110.0
      sm = only(objectsOf(store, :SynchronousMachine))
      @test enumval(sm, :operatingMode) == "SynchronousMachineOperatingMode.generator"
      @test boolval(sm, :earthing) == false
    end

    @testset "rdf:about overlays EQ object (SSH)" begin
      sm = only(objectsOf(store, :SynchronousMachine))
      @test num(sm, :p) == -120.0                    # from SSH overlay
      @test num(sm, :ratedS) == 150.0                # EQ value kept
      @test countOf(store, :SynchronousMachine) == 1 # no duplicate object
    end

    @testset "unresolved references reported" begin
      unresolved = unresolvedReferences(store)
      @test any(u -> u.target == "_missing_container", unresolved)
    end

    @testset "summarizeCGMES on synthetic set" begin
      dir = mktempdir()
      write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
      s = summarizeCGMES(path = dir)
      @test s.object_count == 3
      @test s.unresolved_count >= 1
      @test !s.boundary_missing_hint                 # dangling ref is not a topology class
      @test (:ACLineSegment => 1) in s.class_histogram
    end

    @testset "Result tables stay aligned with long names" begin
      # CGMES bus and branch identifiers routinely exceed the column widths;
      # @sprintf pads but never truncates, so without fitting every following
      # column would shift (reported from a MicroGrid Assembled run).
      @test Sparlectra._fitColumn("short", 25) == "short"
      @test length(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25)) == 25
      @test endswith(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25), "…")
      @test Sparlectra._fitColumn("exactly_twenty_five_chars", 25) == "exactly_twenty_five_chars"
    end

    # The Web UI docs reader serves an allowlist — the CGMES page and the
    # contextual help for its options must be reachable from the interface.
    @testset "Web UI documentation wiring" begin
      page = Sparlectra.resolve_webui_doc_page("cgmes_import")
      @test page !== nothing && page.file == "cgmes_import.md"
      @test isfile(joinpath(dirname(@__DIR__), "docs", "src", page.file))
      # The page must carry the option reference the docs reader links to.
      text = read(joinpath(dirname(@__DIR__), "docs", "src", page.file), String)
      for key in ("cgmes_import.path", "cgmes_import.base_mva", "cgmes_import.require_boundary", "cgmes_import.tap_control", "cgmes_import.ignore_connected")
        @test occursin(key, text)
      end
    end

    # A single ReliCapGrid model is one area of a multi-area system, so its
    # border nodes hang free when imported alone. The combined aliases fetch
    # several areas plus their shared boundary files as ONE delivery. Table
    # consistency is checked here without touching the network; the fetch itself
    # needs GitHub and is exercised by the fixture block below when cached.
    @testset "ReliCapGrid combined aliases" begin
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
      # Aliases must be discoverable through the public alias list.
      all_aliases = Sparlectra.CGMESImporter.allCGMESTestSetAliases()
      @test "relicapgrid_cgm" in all_aliases
      @test "svedala_neighbours" in all_aliases
      @test "portheim" in all_aliases
    end

    # --- real fixtures (optional, offline-safe) -----------------------------
    cache = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(@__DIR__), "data", "CGMES"))
    bc = joinpath(cache, "extracted", "MicroGrid", "BaseCase_BC")
    be = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2")
    asm = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_Assembled_v2")
    bd = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
    if isdir(be) && isdir(bd)
      @testset "MicroGrid BC BE (ENTSO-E fixture)" begin
        s = summarizeCGMES(path = [be, bd])
        @test s.version == "2.4.15"
        @test s.unresolved_count == 0
        @test !s.boundary_missing_hint
        hist = Dict(s.class_histogram)
        @test hist[:TopologicalNode] == 12
        @test hist[:ACLineSegment] == 7
        @test hist[:PowerTransformerEnd] == 9
        # boundary detection: without BD the hint must fire
        s2 = summarizeCGMES(path = be)
        @test s2.unresolved_count > 0
        @test s2.boundary_missing_hint
      end

      # Stage-1 exit criterion: converge and reproduce the SV profile.
      # vm tolerance 2e-4 (the SvVoltage values carry ~6 significant digits,
      # so ~1.5e-4 is the data noise floor), va tolerance 0.05°.
      @testset "Stage 1: import + runpf! + compareWithSV" begin
        res = importCGMES(path = [be, bd], name = "MicroGrid_BE")
        @test length(res.net.nodeVec) == 12
        @test res.slack_bus != ""
        # short-circuit harvest (§7.7): read, not evaluated
        @test length(res.shortcircuit.ac_line_segments) == 7
        @test length(res.shortcircuit.synchronous_machines) == 2
        @test !isempty(res.shortcircuit.equivalent_injections)
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        @test cmp.n == 11
        @test cmp.max_dvm < 2e-4
        @test cmp.max_dva < 0.05
        # SvPowerFlow comparison (D5b): the conformity sets carry injection
        # terminals only; tolerances reflect the SV data noise floor
        @test cmp.flows.n >= 12
        @test cmp.flows.max_dp < 3.0
        @test cmp.flows.max_dq < 2.0
        @test any(r -> r.kind == :shunt, cmp.flows.rows)
        @test any(r -> r.kind == :units, cmp.flows.rows)

        # boundary error path
        @test_throws ErrorException importCGMES(path = be)
      end

      if isdir(asm)
        @testset "Stage 1: assembled model (both sides + X-node EI skip)" begin
          res = importCGMES(path = [asm, bd], name = "MicroGrid_Assembled")
          @test count(m -> occursin("assembled boundary node", m), res.messages) == 10
          @test count(m -> occursin("bus link", m), res.messages) >= 1   # NL busbar coupler
          _, erg = runpf!(res.net, 30, 1e-8, 0)
          @test erg == 0
          cmp = compareWithSV(res)
          @test cmp.max_dvm < 2e-4
          @test cmp.max_dva < 0.05
        end
      end

      # WebUI text-entry path: alias → combined single ZIP (base case plus
      # boundary) in the case directory; importing that ZIP must reproduce
      # the folder import exactly. Offline-safe: the extracted cache exists
      # here, so no download happens.
      @testset "Test-set fetch: alias to combined ZIP" begin
        out = mktempdir()
        z = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = out)
        @test isfile(z) && basename(z) == "cgmes_microgrid_be.zip"
        rz = importCGMES(path = z, name = "from_zip")
        @test length(rz.net.nodeVec) == 12
        _, erg = runpf!(rz.net, 30, 1e-8, 0)
        @test erg == 0
        @test compareWithSV(rz).max_dvm < 2e-4
        @test_throws ErrorException Sparlectra.CGMESImporter.fetchCGMESTestSet("no_such_set"; outdir = out)
        # The WebUI runs through the case selector and the service resolver,
        # not through run_sparlectra_api directly — both must accept the ZIP
        # (regression: the delivery imported fine but was invisible in the
        # selector and rejected by the resolver).
        @test Sparlectra._webui_is_user_selectable_case(z)
        @test basename(Sparlectra._resolve_powerflow_casefile(basename(z), out)) == basename(z)
        @test basename(Sparlectra._resolve_powerflow_casefile(z, out)) == basename(z)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile("does_not_exist.zip", out)

        # Upload-time completeness feedback: the user must learn at import
        # time whether a delivery can be computed, not first at run time.
        sgdir = joinpath(cache, "extracted", "SmallGrid", "BusBranch")
        basezip = joinpath(sgdir, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0.zip")
        bndzip = joinpath(sgdir, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0.zip")
        if isfile(basezip) && isfile(bndzip)
          lonely = mktempdir()
          cp(basezip, joinpath(lonely, "base.zip"))
          # With the local test-set cache present the boundary is supplied
          # automatically; the message must say so and the file must appear.
          lonely_role = Sparlectra._webui_cgmes_upload_role(joinpath(lonely, "base.zip"), lonely)
          @test occursin("boundary supplied automatically", lonely_role) || occursin("BOUNDARY SET MISSING", lonely_role)
          if occursin("boundary supplied automatically", lonely_role)
            @test length(filter(n -> occursin(r"(?i)boundary", n), readdir(lonely))) == 1
            rl = run_sparlectra_api(casefile = joinpath(lonely, "base.zip"), output_dir = mktempdir())
            @test rl.status == :succeeded
          end
          # the pre-analysis must state version, profiles and element counts
          @test occursin("CGMES 2.4.15", lonely_role)
          @test occursin("nodes", lonely_role) && occursin("lines", lonely_role) && occursin("transformers", lonely_role)
          together = mktempdir()
          cp(basezip, joinpath(together, "base.zip"))
          cp(bndzip, joinpath(together, "boundary.zip"))
          @test occursin("ready to compute", Sparlectra._webui_cgmes_upload_role(joinpath(together, "base.zip"), together))
          @test occursin("boundary set", Sparlectra._webui_cgmes_upload_role(joinpath(together, "boundary.zip"), together))
          # and the run picks the neighbouring boundary up automatically
          rr = run_sparlectra_api(casefile = joinpath(together, "base.zip"), output_dir = mktempdir())
          @test rr.status == :succeeded
          @test rr.metadata["cgmes_buses"] == 118
        end

        # Boundary-only deliveries are companions, not runnable cases:
        # they must stay out of the case selector.
        bonly = filter(n -> occursin(r"(?i)boundary", n), readdir(out))
        for b in bonly
          @test !Sparlectra._webui_is_user_selectable_case(joinpath(out, b))
        end
        @test Sparlectra._webui_is_user_selectable_case(z)   # the combined set stays selectable

        # WebUI resolve handler on the cgmes: prefix
        resp = Sparlectra.handle_powerflow_case_resolve(Dict("casefile" => "cgmes:microgrid_be"); output_root = mktempdir(), case_directory = out)
        @test resp.status == 303
        headers = resp.headers isa AbstractDict ? resp.headers : Dict(resp.headers)
        @test occursin("cgmes_microgrid_be.zip", get(headers, "Location", ""))
      end

      # D4b: the framework/API path — CGMES delivery as casefile, boundary
      # set from the cgmes_import config block, dispatch by auto-detection
      @testset "D4b: run_sparlectra_api dispatch on a CGMES delivery" begin
        out = mktempdir()
        cfgfile = joinpath(out, "cfg.yaml")
        write(cfgfile, "cgmes_import:\n  path: \"" * bd * "\"\n")
        r = run_sparlectra_api(casefile = be, config_file = cfgfile, output_dir = out)
        @test r.status == :succeeded
        @test r.metadata["input_format_detected"] == "cgmes"
        @test r.metadata["cgmes_version"] == "2.4.15"
        @test r.metadata["cgmes_buses"] == 12
        @test r.metadata["cgmes_slack_bus"] == "BE-Busbar_4"
        @test r.metadata["cgmes_messages"] > 0
        # a dedicated cgmes.log documents the import next to run.log
        logfile = joinpath(out, "cgmes.log")
        @test isfile(logfile)
        text = read(logfile, String)
        @test occursin("# CGMES import report", text)
        @test occursin("## Class histogram", text)
        @test occursin("## Network built", text)
        @test occursin("## Short-circuit source data", text)
        @test occursin("## Importer messages", text)
      end

      # The report must survive a failed solve — it is written right after the
      # import, not after the power flow (regression: a diagnose run whose
      # power flow failed produced no cgmes.log at all).
      @testset "cgmes.log is written even when the run fails" begin
        out2 = mktempdir()
        cfg2 = joinpath(out2, "c.yaml")
        write(cfg2, "power_flow:\n  max_iter: 1\n  tol: 1.0e-14\n")
        rf = run_sparlectra_api(casefile = be, config_file = cfg2, output_dir = out2)
        @test isfile(joinpath(out2, "cgmes.log"))
        @test occursin("# CGMES import report", read(joinpath(out2, "cgmes.log"), String))
      end

      # Stage 2 (Phase E): start from the SSH tap positions and let the
      # outer-loop controllers find the operating point. Note that the
      # reference SvTapStep positions are NOT reproduced exactly — the CGMES
      # target deadbands of this data set are wide (PST: 35 MW, OLTC: 0.5 kV
      # on 10.5 kV), so a controller legitimately stops one step short. The
      # testable properties are: controllers are created, the loop converges,
      # targets are met within their deadband, and the tap moves toward the
      # SV reference rather than away from it.
      @testset "Stage 2: tap controllers from TapChangerControl" begin
        stage1 = importCGMES(path = [be, bd], name = "s1")
        @test isempty(collect(Sparlectra._tap_controllers(stage1.net)))

        res = importCGMES(path = [be, bd], tap_control = true, name = "s2")
        ctrls = collect(Sparlectra._tap_controllers(res.net))
        @test length(ctrls) == 2
        @test any(c -> c.mode == :branch_active_power && c.control_phase, ctrls)
        @test any(c -> c.mode == :voltage && c.control_ratio, ctrls)

        # the OLTC of this set regulates the slack busbar → disabled by the guard
        @test count(m -> occursin("disabled — target bus is", m), res.messages) == 1
        @test count(c -> !c.enabled, ctrls) == 1

        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
        runres = run_sparlectra(net = res.net, config = cfg)
        @test runres.numerical_converged
        cres = latest_control_result(res.net)
        @test cres.status == :converged

        pst = only(filter(c -> c.mode == :branch_active_power, ctrls))
        @test pst.converged
        @test abs(pst.achieved_p_mw - pst.p_target_mw) <= pst.deadband_p_mw

        # Taps are discrete, so the distance to the reference is measured in
        # STEPS, not in degrees: project the controlled angle back onto the
        # changer's step grid (same evaluation the importer uses, so the
        # asymmetrical PST's non-uniform degree spacing cannot distort the
        # result) and compare integers. The solver moves in a constant
        # `phase_step_deg`, i.e. a linear approximation of that grid — if the
        # approximation ever became too coarse, this assertion is where it
        # would surface.
        pstbr = Sparlectra._find_trafo_branch(res.net, pst.trafo)
        ptc = only(objectsOf(res.store, :PhaseTapChangerAsymmetrical))
        low = Int(round(num(ptc, :lowStep, 0.0)))
        high = Int(round(num(ptc, :highStep, 0.0)))
        nearest_step(angle) = argmin(
          [abs(Sparlectra.CGMESImporter._phaseTapRatioShift_atstep(ptc, s)[2] - angle) for s = low:high],
        ) + low - 1
        step_controlled = nearest_step(pstbr.phase_shift_deg)
        step_sv = 16         # SvTapStep position of BE-TR2_1
        @test abs(step_controlled - step_sv) <= 1
      end
    else
      @info "CGMES MicroGrid fixture not cached — skipping ENTSO-E fixture tests (run examples/experimental/cgmes_fetch_testsets.jl to enable)"
    end

    # FullGrid is ENTSO-E's import/export completeness set: it extends
    # MicroGrid T3 so that every class defined in the CGMES profiles appears
    # at least once (HVDC converters, SVC, AsynchronousMachine, all four PST
    # types). It is explicitly NOT an analytical fixture — no SV accuracy is
    # asserted here. Its job is to prove that the reader walks past unknown
    # classes instead of stumbling, and that they are counted, not dropped.
    fg = joinpath(cache, "extracted", "FullGrid")
    fgbb = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_BB_BE_v2")
    fgnb = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_NB_BE_v4")
    fgbd = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_BD_v1")
    if isdir(fgbb) && isdir(fgnb) && isdir(fgbd)
      @testset "FullGrid: reader completeness across all CGMES classes" begin
        s = summarizeCGMES(path = [fgbb, fgbd])
        @test s.version == "2.4.15"
        @test length(s.class_histogram) > 100      # every profile class present
        @test !s.boundary_missing_hint
        hist = Dict(s.class_histogram)
        # classes far outside the Stage-1 mapping set must still be counted
        for cls in (:VsConverter, :DCLineSegment, :StaticVarCompensator, :AsynchronousMachine, :PhaseTapChangerTabular)
          @test haskey(hist, cls)
        end
        # boundary detection still fires on the richer set
        s2 = summarizeCGMES(path = fgnb)
        @test s2.boundary_missing_hint

        # Node-breaker sets ship the TP profile too, so they read as
        # bus-branch today: the NB variant must yield the same network as
        # the BB variant of the same grid (no node-breaker stage needed).
        rbb = importCGMES(path = [fgbb, fgbd], name = "FullGrid_BB")
        rnb = importCGMES(path = [fgnb, fgbd], name = "FullGrid_NB")
        @test length(rbb.net.nodeVec) == length(rnb.net.nodeVec)
        @test length(rbb.net.branchVec) == length(rnb.net.branchVec)
        @test length(rbb.net.linkVec) == length(rnb.net.linkVec)

        # StaticVarCompensator: mapped as a P = 0 reactive injection
        # (MATPOWER-style); FullGrid's placeholder Ω ratings (±0.99 Ω at
        # 225 kV ↔ ±51 GVAr) trigger the plausibility clamp, and the droop
        # slope is reported as ignored
        @test any(m -> occursin("StaticVarCompensator", m) && occursin("clamped", m), rbb.messages)
        @test any(m -> occursin("StaticVarCompensator", m) && occursin("slope", m), rbb.messages)
        @test any(ps -> something(ps.pVal, -1.0) == 0.0 && something(ps.minQ, 0.0) == -1000.0 && something(ps.maxQ, 0.0) == 1000.0, rbb.net.prosumpsVec)
        @test string(Sparlectra.toProSumptionType("STATICVARCOMPENSATOR")) == "Injection"

        # D6: FullGrid populates the ExternalNetworkInjection short-circuit
        # attributes that MicroGrid lacks — harvested, not evaluated
        eni = only(rbb.shortcircuit.external_network_injections)
        @test eni.maxInitialSymShCCurrent_A !== nothing
        @test eni.minInitialSymShCCurrent_A !== nothing
        @test eni.maxR1ToX1Ratio !== nothing
        @test eni.maxZ0ToZ1Ratio !== nothing
        @test eni.bus !== nothing
      end
    end

    sgbase = joinpath(cache, "extracted", "SmallGrid", "BusBranch")
    sgbc = joinpath(sgbase, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0")
    sgbd = joinpath(sgbase, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0")
    if isdir(sgbc) && isdir(sgbd)
      @testset "Stage 1: SmallGrid bus-branch" begin
        res = importCGMES(path = [sgbc, sgbd], name = "SmallGrid_BB")
        @test length(res.net.nodeVec) == 118
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        @test cmp.max_dvm < 2e-4
        @test cmp.max_dva < 0.05
        @test cmp.flows.n >= 100
        @test cmp.flows.max_dp < 0.5
        @test cmp.flows.max_dq < 0.5
      end
    end

    # --- ReliCapGrid / Svedala: second source, first CGMES 3.0 delivery -----
    #
    # Svedala comes from ENTSO-E's ReliCapGrid repository, not the conformity
    # package, and its layout differs fundamentally: individual CIM/XML files
    # per model, one boundary file per border, plus a shared commonData file
    # (which alone takes the unresolved references from 571 down to 20).
    #
    # Its job in this suite is threefold:
    #   1. it is the only CGMES 3.0 delivery we have — the 3.0 path was
    #      prepared in the schema layer but never exercised by real data;
    #   2. it pins the SSH `Equipment.inService` semantics (a CGMES 3.0
    #      addition): ReliCapGrid parks out-of-service units and switches with
    #      connected=true, inService=false, no SvVoltage — ignoring the flag
    #      imported 413 MW of phantom generation and falsely closed switches
    #      between separately-solved network parts;
    #   3. its SV declares two angle references in two TopologicalIslands, and
    #      with inService honored our electrical-island view AGREES with that
    #      partition — one reference per island under multi_slack.
    #
    # Historical note: this fixture long asserted NO converged power flow
    # (border nodes hanging free). The accumulated fixes changed that — with
    # Equipment.inService honored, per-island references and the single-sided
    # border equivalents standing in for the absent neighbours, Svedala alone
    # now solves AND reproduces its SV state; the testset asserts it.
    svedala = joinpath(cache, "relicapgrid", "cgmes-3.0_ncp-2.5_tc-2.0", "Svedala")
    if isdir(svedala)
      @testset "ReliCapGrid Svedala (CGMES 3.0)" begin
        s = summarizeCGMES(path = svedala)
        @test s.version == "3.0"
        @test !s.boundary_missing_hint            # border files ship with the model
        @test sum(v for (_, v) in s.class_histogram) == 14506
        @test length(s.class_histogram) == 52
        # 20 references stay unresolved even with commonData: they point into
        # the neighbouring areas that this single-area delivery does not carry.
        @test s.unresolved_count == 20

        res = importCGMES(path = svedala, name = "Svedala")
        net = res.net
        # buses are created lazily, so equipment skipped as out of service
        # (32 objects here) never materializes its bus: 132, not the 200 of the
        # raw TP profile. Only 1 of the 47 retained switches is actually in
        # service AND closed.
        @test length(net.nodeVec) == 132
        @test length(net.branchVec) == 150
        @test length(net.linkVec) == 1
        @test count(m -> occursin("out of service (SSH inService=false)", m), res.messages) == 32
        @test any(m -> occursin("VATTENDRAGET_G1", m) && occursin("out of service", m), res.messages)

        # The SV profile names the reference the exporting tool actually used;
        # that outranks referencePriority in the slack chain.
        @test any(m -> occursin("SV declares", m) && occursin("BP_SD-BO3", m) && occursin("HÄLLAN_G1_CNODE", m), res.messages)
        @test any(m -> occursin("slack: taken from the SV angle reference", m), res.messages)

        # The two angle references sit in two different electrical islands —
        # matching the delivery's own two TopologicalIslands. (With inService
        # ignored, falsely-closed switches merged the islands and made the
        # second reference look physically wrong; that was an artifact.)
        bp = net.busDict["BP_SD-BO3"]
        hg = net.busDict["HÄLLAN_G1_CNODE"]
        island_of = Sparlectra.electricalIslandOfBus(net)
        @test island_of[bp] != island_of[hg]
        @test length(Sparlectra.electricalIslandComponents(net)) == 2
        # branch-only view still differs from the electrical view: the one
        # closed link joins two branch-components into one electrical island
        @test length(Sparlectra._active_ac_island_components(net)) == 3

        # multi_slack is the default: each electrical island gets its own
        # SV-declared reference — matching the delivery's TopologicalIslands,
        # and required for the island-wise power flow to run at all (an island
        # without any reference fails _validate_island_references!).
        @test net.slackVec == [bp, hg]
        @test Sparlectra.getNodeType(net.nodeVec[bp]) == Sparlectra.Slack
        @test any(m -> occursin("additional reference", m) && occursin("BP_SD-BO3", m), res.messages)
        # the legacy single-reference behavior stays available
        res_ss = importCGMES(path = svedala, name = "Svedala_single_slack", multi_slack = false)
        @test res_ss.net.slackVec == [hg]
        @test Sparlectra.getNodeType(res_ss.net.nodeVec[res_ss.net.busDict["BP_SD-BO3"]]) != Sparlectra.Slack
        @test any(m -> occursin("BP_SD-BO3", m) && occursin("not promoted to slack", m), res_ss.messages)

        # The single area solves and reproduces its own SV state: the kept
        # single-sided border equivalents stand in for the absent neighbours.
        # (Controlled solver settings — the Q-limit path has its own open
        # solver issues and is not what this fixture is accountable for.)
        ite_sv, erg_sv = runpf!(net, 100, 1e-6, 0; method = :rectangular, islands_enabled = true,
                                islands_diagnostic_continue_after_failure = true,
                                qlimits_enabled = false, opt_flatstart = true)
        @test erg_sv == 0
        cmp_sv = compareWithSV(res)
        @test cmp_sv.n > 100
        @test cmp_sv.rms_dva < 5.0     # measured 2.05°
        @test cmp_sv.rms_dvm < 0.08    # measured 0.037
      end

      # Svedala's five parked generators also carry placeholder voltage targets
      # (RegulatingControl.targetValue = 0.001 kV ≈ 5e-5 pu). With inService
      # honored they are skipped outright — but the plausibility band must
      # still catch such values when the units are forced alive, otherwise a
      # delivery that parks units some OTHER way produces a reference bus at
      # ~0 V and a power flow that "converges" onto an all-zero solution.
      @testset "implausible voltage setpoints are rejected, not obeyed" begin
        # default import: parked units are skipped before their setpoint is read
        res = importCGMES(path = svedala, name = "Svedala_vset")
        @test isempty([m for m in res.messages if startswith(m, "warning:")])

        # ignore_connected revives them — now the band has to fire
        res_ic = importCGMES(path = svedala, name = "Svedala_vset_ic", ignore_connected = true)
        warnings = [m for m in res_ic.messages if startswith(m, "warning:")]
        @test length(warnings) == 5
        @test all(m -> occursin("implausible voltage target", m), warnings)
        # The substituted value has to be derived from the nominal data, and the
        # message must name the unit so the defect is traceable to the source.
        @test any(m -> occursin("VATTENDRAGET_G1", m) && occursin("0.001 kV", m) && occursin("17.0 kV", m), warnings)

        # The affected buses must end up at a usable voltage rather than at zero.
        for name in ("VATTENDRAGET_G1_CNODE", "ATOMSBORG_G1_CNODE", "HÄSTSJÖ_G2_CNODE", "STUPET_G3_CNODE", "RUTHUVUD_G3_CNODE")
          idx = res_ic.net.busDict[name]
          @test res_ic.net.nodeVec[idx]._vm_pu > 0.5
        end

        # The band is a default, not law: widening it accepts the delivery's own
        # values again, which is what a data set that legitimately regulates
        # outside the default band needs.
        res_wide = importCGMES(path = svedala, name = "Svedala_vset_wide", ignore_connected = true, vset_min_pu = 0.0, vset_max_pu = 1.0e9)
        @test isempty([m for m in res_wide.messages if startswith(m, "warning:")])
        @test_throws ArgumentError importCGMES(path = svedala, name = "Svedala_vset_bad", vset_min_pu = 1.2, vset_max_pu = 0.8)
      end

      # The combined seven-area delivery exercises the whole multi-area chain:
      # AC X-node equivalents discarded, the DC border crossing split per side
      # (Stage-0 HVDC), inline VSC converters mapped as fixed PCC injections,
      # and — unlike every single-area model — a power flow that actually
      # converges and reproduces the delivery's own SV state. Offline-safe:
      # fetchCGMESTestSet builds the ZIP from the per-model cache without
      # network when all member folders are present.
      rcg_root = joinpath(cache, "relicapgrid", "cgmes-3.0_ncp-2.5_tc-2.0")
      rcg_members = ("Svedala", "Espheim", "Portheim", "Belgovia", "Galia", "Britheim", "Nordheim")
      if all(isdir(joinpath(rcg_root, m)) for m in rcg_members) && isdir(joinpath(rcg_root, "_shared"))
        @testset "ReliCapGrid combined model (multi-area assembly + Stage-0 HVDC)" begin
          out = mktempdir()
          z = Sparlectra.CGMESImporter.fetchCGMESTestSet("relicapgrid_cgm"; outdir = out)
          res = importCGMES(path = z, name = "ReliCapGrid_CGM", multi_slack = true)
          msgs = res.messages

          # cancelling AC X-node pairs are discarded ...
          @test count(m -> occursin("has both sides present", m), msgs) == 25
          # ... the non-cancelling DC crossing is split per side instead
          @test count(m -> occursin("node split per side", m), msgs) == 1
          @test any(m -> occursin("BP_SD-EH_DC1", m) && occursin("node split per side", m), msgs)
          @test haskey(res.net.busDict, "BP_SD-EH_DC1") && haskey(res.net.busDict, "BP_SD-EH_DC1@2")
          # inline VSC stations become fixed PCC injections
          @test count(m -> occursin("Stage-0 HVDC", m), msgs) == 2

          # the assembled model must solve and reproduce the delivery's SV
          # state. Q-limits stay off here: the active-set path currently cannot
          # handle this case (documented solver follow-up), and the quality
          # gates below are what the import is accountable for.
          net = res.net
          ite, erg = runpf!(net, 100, 1e-6, 0; method = :rectangular, islands_enabled = true,
                            islands_diagnostic_continue_after_failure = true,
                            qlimits_enabled = false, opt_flatstart = true)
          @test erg == 0
          cmp = compareWithSV(res)
          @test cmp.n > 250
          @test cmp.rms_dva < 10.0          # measured 6.3° — angle profile reproduced
          @test cmp.flows.n > 500           # flow comparison works despite isolated buses
          @test cmp.flows.rms_dp < 25.0     # measured 16.3 MW
        end
      end
    else
      @info "ReliCapGrid Svedala not cached — skipping CGMES 3.0 fixture (fetchCGMESTestSet(\"svedala\") enables it)"
    end
  end
end
