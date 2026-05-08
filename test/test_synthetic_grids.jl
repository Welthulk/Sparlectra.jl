# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: test/test_synthetic_grids.jl

function run_synthetic_grid_tests()
  @testset "YAML subset parser" begin
    @test parse_yaml_scalar("true") === true
    @test parse_yaml_scalar("no") === false
    @test parse_yaml_scalar("null") === nothing
    @test parse_yaml_scalar("42") == 42
    @test parse_yaml_scalar("1e-8") ≈ 1e-8
    @test parse_yaml_scalar(":rectangular") === :rectangular
    @test parse_yaml_scalar("'hello # not comment'") == "hello # not comment"
    @test parse_yaml_scalar("[100, 300, 500]") == [100, 300, 500]
    @test as_bool("on") === true
    @test as_bool("off") === false
    @test as_int_vector([4, 9, 16]) == [4, 9, 16]
    @test_throws ArgumentError as_bool("maybe")
    @test_throws ArgumentError as_int_vector([1, 2.5])

    mktemp() do path, io
      write(io, """
      bus_limits: [4, 9, 16]
      solver:
        max_iter: 25
        tol: 1e-8
        method: :rectangular
      synthetic_network:
        aspect_ratio: 1.0 # comment
        label: "demo # grid"
      """)
      close(io)
      cfg = load_yaml_dict(path)
      @test cfg["bus_limits"] == [4, 9, 16]
      @test cfg["solver"]["max_iter"] == 25
      @test cfg["solver"]["tol"] ≈ 1e-8
      @test cfg["solver"]["method"] === :rectangular
      @test cfg["synthetic_network"]["label"] == "demo # grid"
    end

    dst = Dict{String,Any}("solver" => Dict{String,Any}("tol" => 1e-6, "verbose" => 0), "keep" => 1)
    src = Dict{String,Any}("solver" => Dict{String,Any}("tol" => 1e-8), "new" => true)
    merge_yaml_dict!(dst, src)
    @test dst["solver"]["tol"] ≈ 1e-8
    @test dst["solver"]["verbose"] == 0
    @test dst["keep"] == 1
    @test dst["new"] === true

    mktemp() do path, io
      write(io, "root:\n   bad: 1\n")
      close(io)
      @test_throws ErrorException load_yaml_dict(path)
    end
  end

  @testset "synthetic tiled grid builder" begin
    net, meta = build_synthetic_tiled_grid_net(16; aspect_ratio = 1.0, b = 0.01, g = 0.001)
    @test meta.rows >= 2
    @test meta.cols >= 2
    @test meta.actual_buses <= meta.requested_max_buses
    expected = meta.rows * (meta.cols - 1) + (meta.rows - 1) * meta.cols + (meta.rows - 1) * (meta.cols - 1)
    @test meta.branch_count == expected
    @test length(net.nodeVec) == meta.actual_buses
    @test length(net.branchVec) == expected
    @test meta.slack_bus == "B_001_001"
    @test meta.generation_buses == [@sprintf("B_%03d_%03d", meta.rows, 1)]
    @test meta.load_buses == [@sprintf("B_%03d_%03d", 1, meta.cols), @sprintf("B_%03d_%03d", meta.rows, meta.cols)]
    @test synthetic_tiled_grid_bus_index(meta.rows, meta.cols, meta.cols) == meta.actual_buses
    @test geNetBusIdx(net = net, busName = meta.slack_bus) in net.slackVec
    @test any(abs(calcBranchYshunt(branch)) > 0.0 for branch in net.branchVec)

    ybus = createYBUS(net = net, sparse = false, printYBUS = false)
    @test size(ybus) == (meta.actual_buses, meta.actual_buses)
    @test count(!iszero, ybus) > 0

    net2, meta2 = build_synthetic_tiled_grid_net(16; aspect_ratio = 1.0, b = 0.01, g = 0.001)
    @test meta == meta2
    @test net.busDict == net2.busDict

    small_net, small_meta = build_synthetic_tiled_grid_net(9; aspect_ratio = 1.0)
    iterations, erg, _etime = run_net_acpflow(net = small_net, max_ite = 25, tol = 1e-8, verbose = 0, opt_sparse = true, method = :rectangular, show_results = false)
    @test erg == 0
    @test iterations > 0
    @test small_meta.branch_count == small_meta.rows * (small_meta.cols - 1) + (small_meta.rows - 1) * small_meta.cols + (small_meta.rows - 1) * (small_meta.cols - 1)
  end
end
