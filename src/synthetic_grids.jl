# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: src/synthetic_grids.jl

"""
    synthetic_tiled_grid_bus_index(row::Int, col::Int, cols::Int)::Int

Return the deterministic row-major bus index for a synthetic tiled grid.
Rows and columns are one-based, so `row = 1, col = 1` is the upper-left bus.
"""
function synthetic_tiled_grid_bus_index(row::Int, col::Int, cols::Int)::Int
  return (row - 1) * cols + col
end

function _synthetic_tiled_grid_bus_name(row::Int, col::Int)::String
  return @sprintf("B_%03d_%03d", row, col)
end

function _choose_tiled_grid_dims(max_buses::Int, aspect_ratio::Float64)::Tuple{Int,Int}
  max_buses >= 4 || throw(ArgumentError("max_buses must be at least 4 so the tiled grid has at least 2 rows and 2 columns."))
  aspect_ratio > 0.0 || throw(ArgumentError("aspect_ratio must be positive."))

  best_rows = 2
  best_cols = 2
  best_area = 4
  best_score = abs(best_cols / best_rows - aspect_ratio)

  for rows = 2:max_buses
    cols_max = max_buses ÷ rows
    cols_max >= 2 || continue
    for cols = 2:cols_max
      area = rows * cols
      score = abs(cols / rows - aspect_ratio)
      if area > best_area || (area == best_area && score < best_score) || (area == best_area && score == best_score && rows < best_rows)
        best_rows = rows
        best_cols = cols
        best_area = area
        best_score = score
      end
    end
  end

  return best_rows, best_cols
end

function _add_synthetic_pi_line!(net::Net, from_bus::String, to_bus::String, r::Float64, x::Float64, g::Float64, b::Float64)
  @assert from_bus != to_bus "from_bus and to_bus must differ"
  from = geNetBusIdx(net = net, busName = from_bus)
  to = geNetBusIdx(net = net, busName = to_bus)
  vn_kV = getNodeVn(net.nodeVec[from])
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = 1.0, r = r, x = x, b = b, ratedS = nothing, paramsBasedOnLength = false, isPIModel = true)
  acseg.g = g
  push!(net.linesAC, acseg)
  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = 1, values_are_pu = true)
  return nothing
end

"""
    build_synthetic_tiled_grid_net(max_buses::Int; kwargs...) -> (net, metadata)

Build a deterministic synthetic rectangular tiled-grid AC network.

The topology has `rows * cols` buses, horizontal branches, vertical branches,
and one diagonal branch in each rectangular tile from the upper-left bus to the
lower-right bus. Branches are Sparlectra PI-model AC lines whose series
impedance is `r + im*x` in p.u. and whose total shunt admittance is `g + im*b`
in p.u.; Sparlectra's branch model splits that shunt half/half in flow and
Y-bus calculations.

Bus roles follow the other benchmark convention: the upper-left bus is the
slack bus, the lower-left bus has scheduled generation, and the upper-right and
lower-right buses have scheduled PQ loads. Power metadata values are reported
in MW/MVAr; line parameters and voltages are in per unit.

# Keyword defaults
`aspect_ratio = 1.0`, `base_mva = 100.0`, `r = 0.01`, `x = 0.05`,
`g = 0.0`, `b = 0.0`, `load_mw_per_right_corner = 50.0`,
`load_mvar_per_right_corner = 15.0`, `generation_balance = 0.995`,
`vm_slack = 1.0`, `vm_flat = 1.0`.
"""
function build_synthetic_tiled_grid_net(
  max_buses::Int;
  aspect_ratio::Float64 = 1.0,
  base_mva::Float64 = 100.0,
  r::Float64 = 0.01,
  x::Float64 = 0.05,
  g::Float64 = 0.0,
  b::Float64 = 0.0,
  load_mw_per_right_corner::Float64 = 50.0,
  load_mvar_per_right_corner::Float64 = 15.0,
  generation_balance::Float64 = 0.995,
  vm_slack::Float64 = 1.0,
  vm_flat::Float64 = 1.0,
  vn_kV::Float64 = 110.0,
  name::Union{Nothing,String} = nothing,
)
  rows, cols = _choose_tiled_grid_dims(max_buses, aspect_ratio)
  actual_buses = rows * cols
  net_name = isnothing(name) ? "synthetic_tiled_grid_$(rows)x$(cols)" : name
  net = Net(name = net_name, baseMVA = base_mva, flatstart = false)

  for row = 1:rows
    for col = 1:cols
      bus_name = _synthetic_tiled_grid_bus_name(row, col)
      vm = (row == 1 && col == 1) ? vm_slack : vm_flat
      addBus!(net = net, busName = bus_name, vn_kV = vn_kV, vm_pu = vm, va_deg = 0.0, oBusIdx = synthetic_tiled_grid_bus_index(row, col, cols))
    end
  end

  for row = 1:rows
    for col = 1:(cols-1)
      _add_synthetic_pi_line!(net, _synthetic_tiled_grid_bus_name(row, col), _synthetic_tiled_grid_bus_name(row, col + 1), r, x, g, b)
    end
  end
  for row = 1:(rows-1)
    for col = 1:cols
      _add_synthetic_pi_line!(net, _synthetic_tiled_grid_bus_name(row, col), _synthetic_tiled_grid_bus_name(row + 1, col), r, x, g, b)
    end
  end
  for row = 1:(rows-1)
    for col = 1:(cols-1)
      _add_synthetic_pi_line!(net, _synthetic_tiled_grid_bus_name(row, col), _synthetic_tiled_grid_bus_name(row + 1, col + 1), r, x, g, b)
    end
  end

  slack_bus = _synthetic_tiled_grid_bus_name(1, 1)
  generation_bus = _synthetic_tiled_grid_bus_name(rows, 1)
  load_buses = [_synthetic_tiled_grid_bus_name(1, cols), _synthetic_tiled_grid_bus_name(rows, cols)]
  scheduled_load_mw = 2.0 * load_mw_per_right_corner
  scheduled_load_mvar = 2.0 * load_mvar_per_right_corner
  scheduled_generation_mw = generation_balance * scheduled_load_mw
  scheduled_generation_mvar = generation_balance * scheduled_load_mvar

  addProsumer!(net = net, busName = slack_bus, type = "EXTERNALNETWORKINJECTION", vm_pu = vm_slack, va_deg = 0.0, referencePri = slack_bus)
  addProsumer!(net = net, busName = generation_bus, type = "GENERATOR", p = scheduled_generation_mw, q = scheduled_generation_mvar)
  for bus_name in load_buses
    addProsumer!(net = net, busName = bus_name, type = "ENERGYCONSUMER", p = load_mw_per_right_corner, q = load_mvar_per_right_corner)
  end

  result, msg = validate!(net = net)
  result || error("Synthetic tiled-grid network validation failed: $(msg)")

  expected_branch_count = rows * (cols - 1) + (rows - 1) * cols + (rows - 1) * (cols - 1)
  metadata = (
    requested_max_buses = max_buses,
    actual_buses = actual_buses,
    rows = rows,
    cols = cols,
    branch_count = length(net.branchVec),
    expected_branch_count = expected_branch_count,
    slack_bus = slack_bus,
    generation_buses = [generation_bus],
    load_buses = load_buses,
    scheduled_generation_mw = scheduled_generation_mw,
    scheduled_load_mw = scheduled_load_mw,
    scheduled_generation_mvar = scheduled_generation_mvar,
    scheduled_load_mvar = scheduled_load_mvar,
    power_unit = "MW/MVAr",
    line_parameter_unit = "p.u.",
  )
  metadata.branch_count == expected_branch_count || error("Synthetic branch count mismatch: got $(metadata.branch_count), expected $(expected_branch_count)")

  return net, metadata
end

"""
    build_tiled_grid_net(max_buses::Int; kwargs...) -> (net, metadata)

Alias for [`build_synthetic_tiled_grid_net`](@ref).
"""
build_tiled_grid_net(max_buses::Int; kwargs...) = build_synthetic_tiled_grid_net(max_buses; kwargs...)
