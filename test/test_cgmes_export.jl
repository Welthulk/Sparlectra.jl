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
# file: test/test_cgmes_export.jl
# purpose: CGMES export identity and roundtrip contract — structural keys in
#          net.cgmes_ids, parallel-line numbering, deterministic uuid5
#          minting, byte-identical re-export, the duplicate-mRID guard, the
#          power-flow-identical self-roundtrip on EQ+TP+SSH+SV (transformers
#          incl. phase shift and ratio-tap machinery, machines, SVC, loads,
#          shunts, bus links), SSH/SV content, ZIP packaging, tool
#          provenance, and (cache-gated) the MicroGrid roundtrip proving
#          imported mRIDs survive an export, the 3W star reassembly, and the
#          short-circuit evaluation of a re-imported delivery.

using Test
using Sparlectra
using Dates
using UUIDs
using Logging

const _CGMES_EXPORT_CACHE = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(@__DIR__), "data", "CGMES"))

# Fixed header stamp: with `created` pinned the exported files are
# byte-reproducible, which the determinism tests rely on.
const _EXPORT_STAMP = DateTime(2026, 1, 1, 12, 0, 0)

# bare-call defaults are too tight for these fixtures — always solve through
# the explicit signature (maxIte, tol, verbose)
_solve!(net) = runpf!(net, 30, 1e-8, 0)

function _export_test_net()
  net = Net(name = "exportnet", baseMVA = 100.0)
  addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "C", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  # Two parallel lines on the A-B corridor; the second is declared B->A to
  # exercise the lexicographic bus-pair normalization of the structural key.
  addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B", toBus = "A", r_pu = 0.02, x_pu = 0.06, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.015, x_pu = 0.04, b_pu = 0.01, status = 1)
  return net
end

# every stage-2/3 class in one net: two voltage levels, transformer with
# off-nominal ratio + phase shift + ratio-tap machinery, PV machine, slack
# injection, SVC, load, shunt, bus link
function _roundtrip_net()
  net = Net(name = "rt", baseMVA = 100.0)
  addBus!(net = net, busName = "HV1", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "HV1B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "HV2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "MV1", vn_kV = 20.0, vm_pu = 1.0, va_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "HV1", toBus = "HV2", r_pu = 0.01, x_pu = 0.05, b_pu = 0.02, status = 1)
  addLink!(net = net, fromBus = "HV1", toBus = "HV1B", status = 1)
  addPIModelTrafo!(net = net, fromBus = "HV2", toBus = "MV1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, status = 1, ratio = 1.04, shift_deg = 2.0)
  addProsumer!(net = net, busName = "HV1", type = "EXTERNALNETWORKINJECTION", referencePri = "HV1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "HV2", type = "SYNCHRONOUSMACHINE", p = 15.0, q = 0.0, vm_pu = 1.01, isRegulated = true, qMin = -50.0, qMax = 50.0, pMin = 0.0, pMax = 80.0)
  addProsumer!(net = net, busName = "MV1", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
  addShuntMatpower!(net = net, busName = "MV1", Gs = 0.0, Bs = 4.0)
  addProsumer!(net = net, busName = "HV2", type = "STATICVARCOMPENSATOR", p = 0.0, q = 0.0, qMin = -121.0, qMax = 121.0)
  # the concrete class marker the CGMES importer maintains; addProsumer!
  # collapses the type string to the generic Generator component
  net.prosumpsVec[end].comp.cTyp = Sparlectra.StaticVarCompensator
  # ratio-tap machinery on the transformer winding (range export)
  net.trafos[1].side1.taps = Sparlectra.PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 1.1)
  return net
end

function _compare_solved(n1, n2; atol = 1e-6)
  for (bus, i1) in n1.busDict
    i2 = n2.busDict[bus]
    @test isapprox(n1.nodeVec[i1]._vm_pu, n2.nodeVec[i2]._vm_pu; atol = atol)
    @test isapprox(n1.nodeVec[i1]._va_deg, n2.nodeVec[i2]._va_deg; atol = atol)
  end
end

_is_acl_equipment_key(k::AbstractString) = startswith(k, "ACL|") && !endswith(k, "|T1") && !endswith(k, "|T2")

function run_cgmes_export_tests()
  @testset "CGMES export identity" begin
    @testset "structural keys and parallel lines" begin
      net = _export_test_net()
      dir = mktempdir()
      files = writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP)
      @test length(files) == 4
      @test all(isfile, files)
      @test endswith(files[3], "_SSH.xml")
      @test endswith(files[4], "_SV.xml")
      @test haskey(net.cgmes_ids, "MODEL|SSH")
      @test haskey(net.cgmes_ids, "MODEL|SV")
      # tool provenance: version + export date in every file, stamp = created
      for f in files
        head = read(f, String)
        @test occursin("Generated by Sparlectra.jl v$(Sparlectra.version()) on 2026-01-01T12:00:00Z", head)
        @test occursin("<md:Model.description>Sparlectra.jl v$(Sparlectra.version()) export</md:Model.description>", head)
      end
      @test haskey(net.cgmes_ids, "ACL|A|B|1")
      @test haskey(net.cgmes_ids, "ACL|A|B|2")
      @test haskey(net.cgmes_ids, "ACL|B|C|1")
      @test haskey(net.cgmes_ids, "TN|A")
      @test haskey(net.cgmes_ids, "BV|110.0")
      @test haskey(net.cgmes_ids, "MODEL|EQ")
      @test haskey(net.cgmes_ids, "MODEL|TP")
      # parallel lines share the corridor but never an mRID
      @test net.cgmes_ids["ACL|A|B|1"] != net.cgmes_ids["ACL|A|B|2"]
      eq = read(files[1], String)
      tp = read(files[2], String)
      for key in ("ACL|A|B|1", "ACL|A|B|2", "ACL|B|C|1")
        @test occursin("rdf:ID=\"_$(net.cgmes_ids[key])\"", eq)
      end
      @test occursin("rdf:ID=\"_$(net.cgmes_ids["TN|A"])\"", tp)
      # terminal keys follow the equipment sequence (T1 = from side)
      @test haskey(net.cgmes_ids, "ACL|A|B|1|T1")
      @test haskey(net.cgmes_ids, "ACL|A|B|1|T2")
      # the SV profile carries one voltage per bus
      sv = read(files[4], String)
      @test count("<cim:SvVoltage rdf:ID", sv) == 3
    end

    @testset "minted ids are uuid5 over the key" begin
      net = _export_test_net()
      writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      ns = Sparlectra.CGMESImporter.CGMES_UUID_NAMESPACE
      @test net.cgmes_ids["TN|A"] == string(UUIDs.uuid5(ns, "TN|A"))
      @test net.cgmes_ids["ACL|A|B|2"] == string(UUIDs.uuid5(ns, "ACL|A|B|2"))
    end

    @testset "re-export is byte-identical" begin
      net = _export_test_net()
      f1 = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      f2 = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      for i in eachindex(f1)
        @test read(f1[i]) == read(f2[i])
      end
    end

    @testset "independent builds identical, names carry no identity" begin
      n1 = _export_test_net()
      n2 = _export_test_net()
      f1 = writeCGMESFiles(n1; path = mktempdir(), created = _EXPORT_STAMP)
      f2 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
      # both nets start with an empty cgmes_ids dict -> same minted ids
      @test n1.cgmes_ids == n2.cgmes_ids
      for i in eachindex(f1)
        @test read(f1[i]) == read(f2[i])
      end
      # renaming a component changes only its display name, never its mRID
      n2.linesAC[1].comp.cName = "renamed_line"
      f3 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
      @test n2.cgmes_ids == n1.cgmes_ids
      eq3 = read(f3[1], String)
      @test occursin("rdf:ID=\"_$(n1.cgmes_ids["ACL|A|B|1"])\"", eq3)
      @test occursin(">renamed_line<", eq3)
    end

    @testset "self-roundtrip is power-flow-identical" begin
      original = _roundtrip_net()
      @test _solve!(original)[2] == 0
      exported = _roundtrip_net()
      @test _solve!(exported)[2] == 0
      dir = mktempdir()
      notes = String[]
      files = writeCGMESFiles(exported; path = dir, created = _EXPORT_STAMP, notices = notes)
      @test length(files) == 4
      # every class in the fixture is representable — nothing may be dropped
      @test isempty(notes)
      eq = read(files[1], String)
      @test occursin("<cim:Breaker", eq)
      @test occursin("<cim:StaticVarCompensator", eq)
      @test occursin("<cim:RatioTapChanger", eq)
      @test occursin("<cim:PhaseTapChangerLinear", eq)
      res = importCGMES(path = dir, name = "rt_back")
      net2 = res.net
      @test res.slack_bus == "HV1"
      @test length(net2.nodeVec) == length(original.nodeVec)
      @test length(net2.branchVec) == length(original.branchVec)
      @test length(net2.prosumpsVec) == length(original.prosumpsVec)
      @test length(net2.shuntVec) == length(original.shuntVec)
      @test length(net2.linkVec) == length(original.linkVec)
      # the reconstructed transformer carries ratio AND shift exactly
      tbr = [br for br in net2.branchVec if occursin("_2WT_", br.comp.cName)]
      @test length(tbr) == 1
      @test isapprox(tbr[1].ratio, 1.04; atol = 1e-9)
      @test isapprox(tbr[1].angle, 2.0; atol = 1e-9)
      # SVC survives as a regulated unit with its rating-derived Q limits
      svc = [p for p in net2.prosumpsVec if p.comp.cTyp == Sparlectra.StaticVarCompensator]
      @test length(svc) == 1
      @test svc[1].minQ == -121.0
      @test svc[1].maxQ == 121.0
      # the exported SV state starts the re-import at the solution
      its, erg = _solve!(net2)
      @test erg == 0
      @test its <= 2
      _compare_solved(original, net2)
    end

    @testset "SSH and SV profiles carry the operating point" begin
      net = Net(name = "sshnet", baseMVA = 100.0)
      addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A", type = "SYNCHRONOUSMACHINE", p = 25.0, q = 5.0, referencePri = "A", vm_pu = 1.02, isRegulated = true)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 25.0, q = 5.0)
      addShuntMatpower!(net = net, busName = "B", Gs = 0.0, Bs = 4.0)
      @test _solve!(net)[2] == 0
      files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      ssh = read(files[3], String)
      # machine convention: injection-positive storage flips sign in SSH
      @test occursin("<cim:RotatingMachine.p>-25.0</cim:RotatingMachine.p>", ssh)
      @test occursin("<cim:RotatingMachine.q>-5.0</cim:RotatingMachine.q>", ssh)
      @test occursin("<cim:SynchronousMachine.referencePriority>1</cim:SynchronousMachine.referencePriority>", ssh)
      @test occursin("<cim:EnergyConsumer.p>25.0</cim:EnergyConsumer.p>", ssh)
      @test occursin("<cim:EnergyConsumer.q>5.0</cim:EnergyConsumer.q>", ssh)
      # voltage target in kV of the regulated bus
      @test occursin("<cim:RegulatingControl.targetValue>$(1.02 * 110.0)</cim:RegulatingControl.targetValue>", ssh)
      @test occursin("<cim:ShuntCompensator.sections>1</cim:ShuntCompensator.sections>", ssh)
      eq = read(files[1], String)
      @test occursin("RegulatingControlModeKind.voltage", eq)
      @test occursin("<cim:LinearShuntCompensator.bPerSection>$(4.0 / 110.0^2)</cim:LinearShuntCompensator.bPerSection>", eq)
      # SV: the re-import validates cleanly against the exported state
      res = importCGMES(path = dirname(files[1]), name = "sv_back")
      @test isempty(res.no_sv_buses)
      cmp = Sparlectra.CGMESImporter.compareWithSV(res)
      @test maximum(abs(r.dvm) for r in cmp.rows) < 1e-12
      @test maximum(abs(r.dva) for r in cmp.rows) < 1e-12
      @test !isempty(cmp.flows.rows)
      @test maximum(abs(r.dp) for r in cmp.flows.rows) < 1e-9
      @test maximum(abs(r.dq) for r in cmp.flows.rows) < 1e-9
    end

    @testset "zip packaging re-imports directly" begin
      net = _export_test_net()
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP, zip = true)
      @test length(files) == 5
      @test endswith(files[5], "_CGMES.zip")
      res = importCGMES(path = files[5], name = "zip_back")
      @test length(res.net.nodeVec) == 3
      @test length(res.net.linesAC) == 3
    end

    @testset "duplicate mRID aborts before writing" begin
      net = _export_test_net()
      net.cgmes_ids["TN|A"] = "deadbeef-0000-0000-0000-000000000001"
      net.cgmes_ids["TN|B"] = "deadbeef-0000-0000-0000-000000000001"
      dir = mktempdir()
      err = try
        writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP)
        nothing
      catch e
        e
      end
      @test err isa ErrorException
      msg = sprint(showerror, err)
      @test occursin("TN|A", msg)
      @test occursin("TN|B", msg)
      # the guard fires before any file is opened
      @test isempty(readdir(dir))
    end

    bc = joinpath(_CGMES_EXPORT_CACHE, "extracted", "MicroGrid", "BaseCase_BC")
    if isdir(bc)
      @testset "MicroGrid roundtrip preserves imported mRIDs" begin
        res = importCGMES(path = bc, name = "mg_export")
        net = res.net
        acl_keys = sort([k for k in keys(net.cgmes_ids) if _is_acl_equipment_key(k)])
        @test length(acl_keys) == length(net.linesAC)
        @test !isempty(acl_keys)
        # captured ids are canonical (RDF underscore prefix stripped)
        @test all(!startswith(v, "_") for v in values(net.cgmes_ids))
        notes = String[]
        files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP, notices = notes)
        # links and phase shifts export now — the delivery is model-complete
        @test isempty(notes)
        eq = read(files[1], String)
        tp = read(files[2], String)
        for k in acl_keys
          @test occursin("rdf:ID=\"_$(net.cgmes_ids[k])\"", eq)
        end
        # transformer / load / machine / shunt / 3W identity: every recorded
        # equipment id (captured or minted) appears as an rdf:ID in the EQ
        for (kind, nseg) in (("PT|", 4), ("EC|", 3), ("SM|", 3), ("SH|", 3), ("PT3|", 5))
          ks = [k for k in keys(net.cgmes_ids) if startswith(k, kind) && length(split(k, "|")) == nseg]
          @test !isempty(ks)
          @test count(occursin("rdf:ID=\"_$(net.cgmes_ids[k])\"", eq) for k in ks) == length(ks)
        end
        # the reassembled 3W transformer references three ends
        @test count("<cim:TransformerEnd.endNumber>3<", eq) == 1
        # star buses stay internal to the reassembled 3W transformer — their
        # TN is deliberately absent from the export
        tn_keys = sort([k for k in keys(net.cgmes_ids) if startswith(k, "TN|") && !occursin("AUX3WT", uppercase(k))])
        @test !isempty(tn_keys)
        @test occursin("rdf:ID=\"_$(net.cgmes_ids[first(tn_keys)])\"", tp)
        @test !occursin("AUX3WT", tp)
        # a re-import of the export rebuilds the same electrical model —
        # every branch parameter equal, including the PST angle
        back = importCGMES(path = dirname(files[1]), name = "mg_back")
        @test length(back.net.nodeVec) == length(net.nodeVec)
        @test length(back.net.branchVec) == length(net.branchVec)
        @test length(back.net.prosumpsVec) == length(net.prosumpsVec)
        @test length(back.net.shuntVec) == length(net.shuntVec)
        @test length(back.net.linkVec) == length(net.linkVec)
        @test back.slack_bus == res.slack_bus
        busname_of = n -> Dict(v => k for (k, v) in n.busDict)
        nm1, nm2 = busname_of(net), busname_of(back.net)
        brkey = (nm, br, cnt) -> begin
          a, b = nm[br.fromBus], nm[br.toBus]
          p = a <= b ? (a, b) : (b, a)
          k = get(cnt, p, 0) + 1
          cnt[p] = k
          (p[1], p[2], k)
        end
        c1 = Dict{Tuple{String,String},Int}()
        by1 = Dict(brkey(nm1, br, c1) => br for br in net.branchVec)
        c2 = Dict{Tuple{String,String},Int}()
        for br in back.net.branchVec
          o = by1[brkey(nm2, br, c2)]
          @test isapprox(br.r_pu, o.r_pu; atol = 1e-9)
          @test isapprox(br.x_pu, o.x_pu; atol = 1e-9)
          @test isapprox(br.b_pu, o.b_pu; atol = 1e-9)
          @test isapprox(br.g_pu, o.g_pu; atol = 1e-9)
          @test isapprox(br.ratio, o.ratio; atol = 1e-9)
          @test isapprox(br.angle, o.angle; atol = 1e-9)
          @test br.status == o.status
        end
        # renaming an imported line keeps its original mRID on re-export
        target_id = net.cgmes_ids[acl_keys[1]]
        # rename every line: no line may lose its imported id
        for line in net.linesAC
          line.comp.cName = string(line.comp.cName, "_renamed")
        end
        files2 = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
        eq2 = read(files2[1], String)
        @test occursin("rdf:ID=\"_$(target_id)\"", eq2)
        for k in acl_keys
          @test occursin("rdf:ID=\"_$(net.cgmes_ids[k])\"", eq2)
        end
      end

      @testset "short-circuit evaluation survives the roundtrip" begin
        res = importCGMES(path = bc, name = "mg_sc_rt")
        dir = mktempdir()
        writeCGMESFiles(res.net; path = dir, created = _EXPORT_STAMP, sc_line_data = cgmesLineShortCircuitData(res), sc_source = res.shortcircuit)
        back = importCGMES(path = dir, name = "mg_sc_back")
        for case in (:max, :min)
          sc1 = runShortCircuit!(res; case = case)
          sc2 = runShortCircuit!(back; case = case)
          @test length(sc2.rows) == length(sc1.rows)
          r2 = Dict(String(r.bus) => r for r in sc2.rows)
          for r in sc1.rows
            o = r2[String(r.bus)]
            @test r.contains_defaulted_data == o.contains_defaulted_data
            if isfinite(r.ik_kA) && isfinite(o.ik_kA)
              @test isapprox(r.ik_kA, o.ik_kA; atol = 1e-9)
            end
          end
        end
      end

      @testset "harvested zero-sequence line data" begin
        res = importCGMES(path = bc, name = "mg_sc_export")
        scd = cgmesLineShortCircuitData(res)
        # every entry keys an existing line and carries the complete pair
        @test all(1 <= i <= length(res.net.linesAC) for i in keys(scd))
        if !isempty(scd)
          files = writeCGMESFiles(res.net; path = mktempdir(), sc_line_data = scd, created = _EXPORT_STAMP)
          eq = read(files[1], String)
          @test occursin("ACLineSegment.r0", eq)
          @test occursin("ACLineSegment.x0", eq)
        end
      end

      # WebUI checkbox path: a normal power-flow run with export_cgmes writes
      # the profile artifacts (imported mRIDs preserved, zero-sequence data
      # riding along, plus the re-importable delivery zip) next to the run
      # artifacts.
      @testset "service run with export_cgmes" begin
        root = mktempdir()
        cfg = joinpath(root, "c.yaml")
        write(cfg, "power_flow:\n  max_iter: 40\n")
        # the service path accepts case FILES; use the combined alias ZIP
        # (base + boundary in one delivery, packed from the local cache)
        zipcase = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = mktempdir())
        resp = start_powerflow_run(Dict("casefile" => zipcase, "config_file" => cfg, "output_root" => root, "export_cgmes" => true))
        @test resp["success"] === true
        @test resp["metadata"]["cgmes_export_status"] == "completed"
        rundir = joinpath(root, resp["run_id"])
        # only the combined delivery zip is an artifact — no loose profile files
        @test isempty(filter(f -> endswith(f, ".xml"), readdir(rundir)))
        zips = filter(f -> endswith(f, "_CGMES.zip"), readdir(rundir))
        @test length(zips) == 1
        @test resp["metadata"]["cgmes_export_files"] == only(zips)
        # the delivery zip re-imports directly and must reuse the imported
        # mRIDs and carry the harvested zero-sequence data
        back = importCGMES(path = joinpath(rundir, only(zips)), name = "mg_service_back")
        res = importCGMES(path = zipcase, name = "mg_service_ref")
        acl_keys_ref = sort([k for k in keys(res.net.cgmes_ids) if _is_acl_equipment_key(k)])
        @test !isempty(acl_keys_ref)
        for k in acl_keys_ref
          @test get(back.net.cgmes_ids, k, nothing) == res.net.cgmes_ids[k]
        end
        sc_ref = cgmesLineShortCircuitData(res)
        @test resp["metadata"]["cgmes_export_sc_lines"] == length(sc_ref)
        isempty(sc_ref) || @test !isempty(cgmesLineShortCircuitData(back))
        @test occursin("CGMES export:", read(joinpath(rundir, "run.log"), String))
      end
    else
      @info "CGMES export roundtrip: MicroGrid fixture not cached — skipping (run examples/experimental/cgmes_fetch_testsets.jl to enable)"
    end
  end
end
