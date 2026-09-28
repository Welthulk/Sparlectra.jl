# Copyright 2023–2026 Udo Schmitz
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

# file: test/extended/test_dtf_importer.jl
# purpose: extended tests for the native DTF importer: focused synthetic
#          DTFCase checks, the configured tap-changer model, the typed
#          phase-tap model on windings, and the untracked FOR001 fixture
# The reference network is untracked; both documented drop-in locations are
# accepted (test/fixtures/dtf/ takes precedence, data/DTF/ is where the
# validation examples already read it from).
const _DTF_FIXTURE_CANDIDATES = (joinpath(@__DIR__, "..", "fixtures", "dtf", "FOR001.DAT"), joinpath(@__DIR__, "..", "..", "data", "DTF", "FOR001.DAT"))
const DTF_FIXTURE = let i = findfirst(isfile, _DTF_FIXTURE_CANDIDATES)
  _DTF_FIXTURE_CANDIDATES[i === nothing ? 1 : i]
end

function _synthetic_dtf_case(; kind::Char = 'T', g_s::Float64 = 4.0e-5, b_s::Float64 = -1.0e-5, path::String = "synthetic")
  return Sparlectra.DTFImporter.DTFCase(
    path,
    100.0,
    Sparlectra.DTFImporter.DTFParams("", Float64[]),
    ["synthetic"],
    [110.0],
    Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 0, "SLACK"),
    [Sparlectra.DTFImporter.DTFBranch("", 1, kind, 1, "A", "PV", "SLACK", 1.21, 6.05, g_s, b_s, nothing)],
    Sparlectra.DTFImporter.DTFCompensation[],
    Sparlectra.DTFImporter.DTFTransformerControl[],
    [
      Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
      Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
    ],
    [Sparlectra.DTFImporter.DTFOutage("", 1, 'T', 1, "A", "PV", "SLACK")],
    [Sparlectra.DTFImporter.DTFTrailingRecord("trailing branch echo", 1, :branch, :echo, nothing)],
  )
end

function _synthetic_dtf_case_with_tap_control(; longitudinal_range_percent::Float64, actual_tap_step::Int, max_tap_step::Int, added_voltage_angle_deg::Float64 = 0.0)
  control = Sparlectra.DTFImporter.DTFTransformerControl(
    "", 1, "PV", "SLACK", "A", "", "PV", "SLACK",
    110.0, 110.0, longitudinal_range_percent, added_voltage_angle_deg, max_tap_step, actual_tap_step,
    nothing, nothing, nothing,
  )
  return Sparlectra.DTFImporter.DTFCase(
    "synthetic_tap",
    100.0,
    Sparlectra.DTFImporter.DTFParams("", Float64[]),
    ["synthetic_tap"],
    [110.0],
    Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 1, "SLACK"),
    [Sparlectra.DTFImporter.DTFBranch("", 1, 'T', 1, "A", "PV", "SLACK", 1.21, 6.05, 0.0, 0.0, nothing)],
    Sparlectra.DTFImporter.DTFCompensation[],
    [control],
    [
      Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
      Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
    ],
    Sparlectra.DTFImporter.DTFOutage[],
    Sparlectra.DTFImporter.DTFTrailingRecord[],
  )
end

function run_dtf_importer_tests()
  @testset "native DTF importer focused synthetic checks" begin (function ()
    case = _synthetic_dtf_case()
    branch = only(case.branches)
    pu = Sparlectra.DTFImporter._branch_pu(case, branch)
    zbase = 110.0^2 / 100.0
    @test pu.r ≈ 1.21 / zbase
    @test pu.x ≈ 6.05 / zbase
    @test pu.g ≈ 4.0e-5 * zbase
    @test pu.b ≈ -1.0e-5 * zbase

    net = Sparlectra.DTFImporter.build_net(case)
    @test net isa Sparlectra.Net
    @test length(net.trafos) == 1
    @test length(net.branchVec) == 1
    winding = only(net.trafos).side1
    @test winding.g ≈ pu.g
    @test winding.b ≈ pu.b
    @test net.branchVec[1].g_pu ≈ pu.g
    @test net.branchVec[1].b_pu ≈ pu.b
    @test Sparlectra.getTrafoRXBG(winding) == (winding.r, winding.x, winding.b, winding.g)
    r_pu, x_pu, b_pu, g_pu = Sparlectra.getTrafoRXBG_pu(winding, 110.0, 100.0)
    @test r_pu ≈ pu.r
    @test x_pu ≈ pu.x
    @test b_pu ≈ pu.b
    @test g_pu ≈ pu.g
    @test Sparlectra.bus_shunt_totals_pu(net).total_g_pu ≈ 0.0
    @test get(net.matpower_branch_metadata, 1, nothing).transformer_loss_allocation == :native_branch_pi

    zero_g_net = Sparlectra.DTFImporter.build_net(_synthetic_dtf_case(g_s = 0.0))
    @test zero_g_net.branchVec[1].g_pu == 0.0
    @test Sparlectra.bus_shunt_totals_pu(zero_g_net).total_g_pu ≈ 0.0

    @test length(case.outages) == 1
    @test only(case.outages).from == "PV"
    @test length(case.trailing_records) == 1

    # The classic result print shows `net.name` as "Case"; it must reflect the
    # originating file name, not a free-text description card from the file.
    named_net = Sparlectra.DTFImporter.build_net(_synthetic_dtf_case(path = "/some/dir/FOR001B.DAT"))
    @test named_net.name == "FOR001B.DAT"
    @test net.name == "synthetic"
  end)() end

  # SCF export of a DTF-sourced net (issue #342): the writer is net-based,
  # not format-based, so a DTF import must reach the case format with its
  # reference names and its tap data intact.
  @testset "SCF export of a DTF-sourced net" begin (function ()
    dtf_net = Sparlectra.DTFImporter.build_net(_synthetic_dtf_case_with_tap_control(longitudinal_range_percent = 10.0, actual_tap_step = 7, max_tap_step = 10))
    root = Sparlectra.net_to_scf(dtf_net; source_format = "dtf", source_reference = "synthetic")
    @test root["sparlectra"]["meta"]["source_format"] == "dtf"
    # bus reference names (the DTF station names) survive the export
    names = Set(String(e["name"]) for e in values(root["sparlectra"]["extra"]))
    @test all(bus -> bus in names, keys(dtf_net.busDict))
    # the transformer reaches data.generic_branch with a live ratio
    @test haskey(root["data"], "generic_branch")
    @test all(r -> r["k"] > 0.0, root["data"]["generic_branch"])
    d = mktempdir()
    a = exportSCF(dtf_net; file = joinpath(d, "a.scf.json"), source_format = "dtf")
    b = exportSCF(dtf_net; file = joinpath(d, "b.scf.json"), source_format = "dtf")
    @test read(a, String) == read(b, String)
  end)() end

  @testset "native DTF importer honors configured tap-changer model" begin (function ()
    case = _synthetic_dtf_case_with_tap_control(longitudinal_range_percent = 10.0, actual_tap_step = 7, max_tap_step = 10)

    net_ideal = Sparlectra.DTFImporter.build_net(case; tap_changer_model = :ideal)
    net_corrected = Sparlectra.DTFImporter.build_net(case; tap_changer_model = :impedance_correction)

    branch_ideal = only(net_ideal.branchVec)
    branch_corrected = only(net_corrected.branchVec)

    tap_fraction = (10.0 / 100.0) * 7 / 10
    expected_factor = (1.0 + tap_fraction)^2
    @test isapprox(branch_corrected.r_pu, branch_ideal.r_pu * expected_factor; atol = 1e-12)
    @test isapprox(branch_corrected.x_pu, branch_ideal.x_pu * expected_factor; atol = 1e-12)
    @test get(net_corrected.matpower_branch_metadata, 1, nothing).tap_changer_model == :impedance_correction
    @test isapprox(get(net_corrected.matpower_branch_metadata, 1, nothing).tap_impedance_correction_factor, expected_factor; atol = 1e-12)
    @test get(net_ideal.matpower_branch_metadata, 1, nothing).tap_impedance_correction_factor == 1.0
  end)() end

  # The shipped demo deck: the text path of the importer on a file every
  # checkout has (the decks of the testsets below are local files). Cards
  # in fixed columns, a transformer control with its typed model, an outage
  # record that resolves to one branch.
  @testset "shipped demo deck sp_dtf5" begin (function ()
    deck = joinpath(dirname(@__DIR__), "..", "data", "dtf_demo", "sp_dtf5.DAT")
    case = Sparlectra.DTFImporter.read_dtf(deck)
    @test (length(case.buses), length(case.branches), length(case.transformer_controls), length(case.outages)) == (5, 6, 1, 2)
    @test case.size.slack == "NORD"
    @test case.nominal_voltages_kv == [110.0, 20.0]
    control = only(case.transformer_controls)
    @test (control.from, control.to, control.longitudinal_range_percent, control.max_tap_step, control.actual_tap_step) == ("SUED", "STADT", 16.0, 8, 2)
    net = Sparlectra.DTFImporter.build_net(case)
    @test only(net.trafos).side1.phase_taps !== nothing
    ite, erg = runpf!(net, 50, 1e-8, 0)
    @test erg == 0
    vm = Dict(name => net.nodeVec[idx]._vm_pu for (name, idx) in net.busDict)
    @test isapprox(vm["NORD"], 1.02; atol = 1e-9)
    @test isapprox(vm["SUED"], 1.0050214287; atol = 1e-6)
    @test isapprox(vm["STADT"], 1.0240776060; atol = 1e-6)
    for outage in case.outages
      matches = Sparlectra.DTFImporter.find_outage_branch_indices(case, outage)
      @test length(matches) == 1
      out = Sparlectra.DTFImporter.build_net(case)
      Sparlectra.DTFImporter.apply_single_branch_outage!(out, only(matches))
      _, out_erg = runpf!(out, 50, 1e-8, 0)
      @test out_erg == 0
    end
    # the front doors read the same deck, and no format has to be named
    cfg = load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
    imported = Sparlectra.import_case(deck, cfg)
    @test imported.format === :dtf_for001
    @test length(imported.net.nodeVec) == 5 && length(imported.net.branchVec) == 6
    run = redirect_stdout(devnull) do
      run_sparlectra(casefile = deck)
    end
    @test run isa Sparlectra.SparlectraRunResult
  end)() end

  # A DTF deck is recognised by its content, and the reader is the judge:
  # what `read_dtf` takes as a network is a deck under any name, everything
  # else is refused. The texts below are written here; the last rows are
  # the local reference files, which are in no repository.
  @testset "format detection by content" begin (function ()
    deck_text = read(joinpath(dirname(@__DIR__), "..", "data", "dtf_demo", "sp_dtf5.DAT"), String)
    cards = split(deck_text, '\n')
    d = mktempdir()
    file = (name, text) -> (path = joinpath(d, name); write(path, text); path)
    # the same deck under a name that says nothing
    @test Sparlectra.DTFImporter.is_dtf_deck(file("network.txt", deck_text))
    @test Sparlectra._detect_case_format(joinpath(d, "network.txt")) === :dtf_for001
    refused = [
      ("plain text", "plain unsupported data\n"),
      ("numbers", "1 2 3\n4 5 6\n"),
      ("a result report", "LASTFLUSSERGEBNIS FOR002\nKNOTEN   U/KV   WINKEL     P/MW    Q/MVAR\nNORD    112.200    0.000   48.484   12.100\nWEST    111.132   -1.015  -30.000  -10.000\nOST     112.000   -0.406   35.000    4.300\nSUED    110.552   -1.529  -25.000   -8.000\nSTADT    20.482   -4.381  -18.000   -6.000\nVERLUSTE 0.484 MW\n"),
      ("an outage list without a network", "AUSFALL\nL1   WEST      OST     \nENDE\n"),
      # the size card counts a bus more than the deck has
      ("a deck shorter than its size card says", join(vcat(cards[1:6], ["    6    6    0    1  NORD"], cards[8:end]), '\n')),
      # no bus card of type 2: the network has no reference
      ("a deck without a slack bus card", replace(deck_text, "21   NORD" => "01   NORD")),
    ]
    for (i, (label, text)) in enumerate(refused)
      path = file("refused_$(i).DAT", text)
      @testset "$(label)" begin
        @test !Sparlectra.DTFImporter.is_dtf_deck(path)
        @test_throws ArgumentError Sparlectra._detect_case_format(path)
      end
    end
    # a named format is still taken as named
    @test Sparlectra._detect_case_format(joinpath(d, "refused_1.DAT"); requested = :dtf_for001) === :dtf_for001
    local_dir = joinpath(dirname(@__DIR__), "..", "data", "DTF")
    for (name, expected) in [("FOR001.DAT", true), ("FOR001B.DAT", true), ("FOR001C.DAT", true), ("FOR001D.DAT", true), ("FOR001E.DAT", true), ("FOR002.DAT", false)]
      path = joinpath(local_dir, name)
      if isfile(path)
        println("      format detection: ", name, " RAN")
        @test Sparlectra.DTFImporter.is_dtf_deck(path) == expected
      else
        println("      format detection: ", name, " SKIPPED (not present under data/DTF)")
      end
    end
  end)() end

  @testset "native DTF full local fixture" begin
    if !isfile(DTF_FIXTURE)
      @info "Skipping full FOR001 fixture validation; external reference network is not tracked. Place a local file at test/fixtures/dtf/FOR001.DAT or data/DTF/FOR001.DAT for manual validation."
      @test_skip "FOR001 full reference fixture not available"
      return
    end
    case = Sparlectra.DTFImporter.read_dtf(DTF_FIXTURE)
    @test case.size.NGES == 13
    @test length(case.branches) == 27
    @test count(b -> b.kind == 'T', case.branches) == 5
    net = Sparlectra.createNetFromDTFFile(DTF_FIXTURE)
    @test length(net.branchVec) == 27
    @test Sparlectra.bus_shunt_totals_pu(net).total_g_pu ≈ 0.0
  end

  @testset "DTF persists the typed phase-tap model on the winding" begin
    # X(α)-coupling precondition: the transiently built PhaseTapChangerModel
    # must survive onto winding.phase_taps for controlled transformers
    # (FOR001E carries skew/longitudinal tap controls). Pure attachment —
    # the numeric branch results are guarded by the DTF validation suites.
    fixture = joinpath(dirname(@__DIR__), "..", "data", "DTF", "FOR001E.DAT")
    if !isfile(fixture)
      @info "FOR001E.DAT not available — skipping phase-tap persistence check"
      @test_skip "FOR001E fixture not available"
      return
    end
    net = Sparlectra.createNetFromDTFFile(fixture)
    models = [tf.side1.phase_taps for tf in net.trafos if tf.side1.phase_taps !== nothing]
    @test !isempty(models)
    @test all(m -> m.kind === :asymmetrical, models)
    @test all(m -> m.winding_connection_angle_deg !== nothing, models)
    # the skew transformer carries its regulating-vector connection angle
    @test any(m -> m.winding_connection_angle_deg != 0.0, models)
  end
end
