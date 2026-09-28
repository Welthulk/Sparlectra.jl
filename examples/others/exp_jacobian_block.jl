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

# Date: 2026-09-28
# file: examples/others/exp_jacobian_block.jl
# purpose: reads the 2x2 block of the rectangular Newton-Raphson Jacobian
#          that differentiates the equations of one bus by the voltage
#          variables of another, and shows where the block comes from:
#          a diagram of the network (SVG and console) with the bus types
#          and the position k of every bus in the Jacobian, a check of the
#          block against finite differences of the mismatch, and the four
#          entries over the Newton iterations as a table

using Sparlectra
using LinearAlgebra
using SparseArrays
using Printf

include(joinpath(@__DIR__, "example_header.jl"))

const OUT_DIR = joinpath(dirname(@__DIR__), "_out", "jacobian_block")

"""
    jacobianBlock(J, net, eqBus, varBus; slack_idx = net.slackVec[1])

The 2x2 block of the rectangular Jacobian `J` that holds the derivatives of
the two equations of bus `eqBus` by the two voltage variables of bus
`varBus`, as `(rows, cols, block)`.

Layout of `J` (the one of `mismatch_rectangular` and
`build_rectangular_jacobian_pq_pv_sparse`), with `k` the position of a bus
in the list of the non-slack buses in network order:

- rows, two per non-slack bus, interleaved: `2k-1` is the active-power
  equation, `2k` the reactive-power equation of a PQ bus or the voltage
  equation of a PV bus;
- columns: `k` is `Vr_k`, `(n-1)+k` is `Vi_k`.

With a distributed slack `J` has one more row and column at its end, which
leaves this layout untouched. The block therefore reads

    [ dP/dVr      dP/dVi
      d(Q|V)/dVr  d(Q|V)/dVi ]

`J` has to belong to the whole network with one reference bus: a network
that solves in several islands, or with buses merged by links, has a
Jacobian per island on its own numbering. The slack bus has no rows and
columns and is refused.
"""
function jacobianBlock(J::AbstractMatrix, net::Net, eqBus::String, varBus::String; slack_idx::Int = net.slackVec[1])
  n = length(net.nodeVec)
  size(J, 1) in (2 * (n - 1), 2 * (n - 1) + 1) || error("the dimension of J ($(size(J, 1))) does not fit a network of $(n) buses")

  i = geNetBusIdx(net = net, busName = eqBus)
  j = geNetBusIdx(net = net, busName = varBus)
  (i == slack_idx || j == slack_idx) && error("the slack bus has no rows and columns in J")

  # position in the list of the non-slack buses
  pos(b) = b < slack_idx ? b : b - 1

  rP = 2 * pos(i) - 1
  rQV = rP + 1
  cVr = pos(j)
  cVi = (n - 1) + pos(j)

  return (rows = (rP, rQV), cols = (cVr, cVi), block = Matrix(J[[rP, rQV], [cVr, cVi]]))
end

# The inputs of the mismatch and the Jacobian at the state the network holds.
function jacobian_inputs(net::Net)
  Ybus = sparse(createYBUS(net = net))
  V = ComplexF64[something(node._vm_pu, 1.0) * cis(deg2rad(something(node._va_deg, 0.0))) for node in net.nodeVec]
  bus_types, Vset, slack_idx = Sparlectra.extract_bus_types_and_vset(net)
  S = Sparlectra.buildComplexSVec(net)
  return (; Ybus, V, S, bus_types, Vset, slack_idx)
end

jacobian_at(p, V) = Sparlectra.build_rectangular_jacobian_pq_pv_sparse(p.Ybus, V, p.bus_types, p.Vset, p.slack_idx)
mismatch_at(p, V) = Sparlectra.mismatch_rectangular(p.Ybus, V, p.S, p.bus_types, p.Vset, p.slack_idx)

# The block by central differences of the mismatch: column by column the
# real and the imaginary part of the voltage of varBus are moved by h.
function finite_difference_block(p, net::Net, eqBus::String, varBus::String; h::Float64 = 1e-6)
  ref = jacobianBlock(jacobian_at(p, p.V), net, eqBus, varBus; slack_idx = p.slack_idx)
  j = geNetBusIdx(net = net, busName = varBus)
  block = zeros(2, 2)
  for (c, step) in enumerate((h + 0.0im, h * im))
    up = copy(p.V)
    down = copy(p.V)
    up[j] += step
    down[j] -= step
    dF = (mismatch_at(p, up) .- mismatch_at(p, down)) ./ (2h)
    block[:, c] = dF[collect(ref.rows)]
  end
  return block
end

# Newton iterations from a flat start, with the block recorded at every
# iterate. The bus types stay as the solved network has them (no limit
# switching), so the layout of J is the same in every step.
function block_history(p, net::Net, eqBus::String, varBus::String; iterations::Int = 6)
  n = length(p.V)
  V = ComplexF64[p.bus_types[k] == :PQ ? 1.0 + 0.0im : p.Vset[k] + 0.0im for k in 1:n]
  V[p.slack_idx] = p.V[p.slack_idx]
  non_slack = [k for k in 1:n if k != p.slack_idx]
  history = Vector{NamedTuple}()
  for it in 0:iterations
    J = jacobian_at(p, V)
    F = mismatch_at(p, V)
    push!(history, (iteration = it, mismatch = maximum(abs, F), block = jacobianBlock(J, net, eqBus, varBus; slack_idx = p.slack_idx).block))
    step = -(J \ F)
    for (k, bus) in enumerate(non_slack)
      V[bus] += step[k] + im * step[(n - 1) + k]
    end
  end
  return history
end

# The network as the Jacobian sees it: every bus with its type and its
# position k in the list of the non-slack buses, every connection between
# two buses with the number of parallel systems, a transformer where the two
# buses differ in nominal voltage.
function network_structure(net::Net, slack_idx::Int)
  n = length(net.nodeVec)
  names = Dict{Int,String}(idx => bus for (bus, idx) in net.busDict)
  pos(b) = b < slack_idx ? b : b - 1
  buses = [(index = i, name = names[i], vn = Sparlectra.getNodeVn(net.nodeVec[i]), type = i == slack_idx ? "Slack" : string(Sparlectra.getNodeType(net.nodeVec[i])), k = i == slack_idx ? 0 : pos(i)) for i in 1:n]
  systems = Dict{Tuple{Int,Int},Int}()
  for br in net.branchVec
    br.status == 1 || continue
    pair = minmax(Int(br.fromBus), Int(br.toBus))
    systems[pair] = get(systems, pair, 0) + 1
  end
  edges = [(from = a, to = b, systems = count, transformer = buses[a].vn != buses[b].vn) for ((a, b), count) in sort!(collect(systems); by = first)]
  return (buses = buses, edges = edges)
end

# The structure on the console: one line per bus with its neighbours.
# "==" is a connection of two systems, "-T-" a transformer, "--" a line.
function print_network(structure, marked::Vector{String})
  width = maximum(length(b.name) for b in structure.buses)
  println("network (k = position of the bus in the Jacobian, * = a bus of the example):")
  for b in structure.buses
    links = String[]
    for e in structure.edges
      (e.from == b.index || e.to == b.index) || continue
      other = structure.buses[e.from == b.index ? e.to : e.from]
      push!(links, string(e.transformer ? "-T- " : e.systems > 1 ? "== " : "-- ", other.name))
    end
    label = b.k == 0 ? "Slack      " : @sprintf("%-5s k=%-2d ", b.type, b.k)
    println("  ", b.name in marked ? "*" : " ", " ", rpad(b.name, width), "  ", label, " ", join(links, "   "))
  end
  println("  rows of J: 2k-1 (dP), 2k (dQ or dV); columns: k (Vr), ", length(structure.buses) - 1, "+k (Vi)")
  return nothing
end

# Placement for the diagram: the buses in layers by their distance from the
# slack bus (so a connection joins neighbouring layers or two buses of one
# layer), ordered inside a layer by the mean position of their neighbours
# in the layer before, which keeps the crossings few.
function diagram_layout(structure, slack_idx::Int)
  n = length(structure.buses)
  neighbours = [Int[] for _ in 1:n]
  for e in structure.edges
    push!(neighbours[e.from], e.to)
    push!(neighbours[e.to], e.from)
  end
  depth = fill(-1, n)
  depth[slack_idx] = 0
  queue = [slack_idx]
  while !isempty(queue)
    b = popfirst!(queue)
    for o in neighbours[b]
      depth[o] == -1 || continue
      depth[o] = depth[b] + 1
      push!(queue, o)
    end
  end
  last_layer = maximum(depth)
  for i in 1:n
    depth[i] == -1 && (depth[i] = last_layer + 1)   # not connected to the slack bus
  end
  row = zeros(Float64, n)
  for d in 0:maximum(depth)
    layer = [i for i in 1:n if depth[i] == d]
    key = i -> begin
      before = [row[o] for o in neighbours[i] if depth[o] == d - 1]
      isempty(before) ? Float64(i) : sum(before) / length(before)
    end
    sort!(layer; by = i -> (key(i), i))
    for (r, i) in enumerate(layer)
      row[i] = r
    end
  end
  return (depth = depth, row = row)
end

# The diagram as plain SVG: a box per bus (name, type, k), a line per
# connection, two lines for two systems, two circles for a transformer.
function write_network_diagram(path::String, structure, slack_idx::Int, marked::Vector{String})
  layout = diagram_layout(structure, slack_idx)
  bw, bh, dx, dy, margin = 150, 46, 215, 92, 30
  rows_of(d) = count(==(d), layout.depth)
  tallest = maximum(rows_of(d) for d in 0:maximum(layout.depth))
  # every layer is centred on the tallest one
  cx(i) = margin + bw / 2 + dx * layout.depth[i]
  cy(i) = margin + bh / 2 + dy * (layout.row[i] - 1 + (tallest - rows_of(layout.depth[i])) / 2)
  W = round(Int, 2 * margin + bw + dx * maximum(layout.depth))
  H = round(Int, 2 * margin + bh + dy * (tallest - 1) + 84)
  fill_of(b) = b.type == "Slack" ? ("#6e2a14", "#f3b9a5") : b.type == "PV" ? ("#0b5345", "#9de0cf") : ("#454540", "#dcdcd4")
  r1(v) = round(v; digits = 1)
  mkpath(dirname(path))
  open(path, "w") do io
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(W)\" height=\"$(H)\" font-family=\"sans-serif\">")
    println(io, "<rect width=\"$(W)\" height=\"$(H)\" fill=\"white\"/>")
    for e in structure.edges
      x1, y1, x2, y2 = cx(e.from), cy(e.from), cx(e.to), cy(e.to)
      len = hypot(x2 - x1, y2 - y1)
      ox, oy = -(y2 - y1) / len, (x2 - x1) / len      # unit normal, for the second system
      same_layer = layout.depth[e.from] == layout.depth[e.to]
      far = same_layer && abs(layout.row[e.from] - layout.row[e.to]) > 1
      for shift in (e.systems > 1 ? (-4.0, 4.0) : (0.0,))
        if far
          # two buses of one layer with others between them: an arc beside the boxes
          bend = bw / 2 + 26 + 10 * abs(layout.row[e.from] - layout.row[e.to])
          println(io, "<path d=\"M $(r1(x1 + bw / 2)) $(r1(y1 + shift)) Q $(r1(x1 + bend + shift)) $(r1((y1 + y2) / 2)) $(r1(x2 + bw / 2)) $(r1(y2 - shift))\" fill=\"none\" stroke=\"#8a8a80\" stroke-width=\"1.6\"/>")
        else
          println(io, "<line x1=\"$(r1(x1 + shift * ox))\" y1=\"$(r1(y1 + shift * oy))\" x2=\"$(r1(x2 + shift * ox))\" y2=\"$(r1(y2 + shift * oy))\" stroke=\"#8a8a80\" stroke-width=\"1.6\"/>")
        end
      end
      if e.transformer
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        far && (mx = x1 + bw / 2 + 0.5 * (bw / 2 + 26 + 10 * abs(layout.row[e.from] - layout.row[e.to]) - bw / 2) + 6)
        ux, uy = (x2 - x1) / len, (y2 - y1) / len
        for side in (-7.0, 7.0)
          println(io, "<circle cx=\"$(r1(mx + side * ux))\" cy=\"$(r1(my + side * uy))\" r=\"9\" fill=\"white\" stroke=\"#8a8a80\" stroke-width=\"1.6\"/>")
        end
      end
    end
    for b in structure.buses
      x, y = cx(b.index) - bw / 2, cy(b.index) - bh / 2
      background, text = fill_of(b)
      frame = b.name in marked ? " stroke=\"#e8a100\" stroke-width=\"3\"" : ""
      println(io, "<rect x=\"$(r1(x))\" y=\"$(r1(y))\" width=\"$(bw)\" height=\"$(bh)\" rx=\"6\" fill=\"$(background)\"$(frame)/>")
      println(io, "<text x=\"$(r1(cx(b.index)))\" y=\"$(r1(y + 19))\" text-anchor=\"middle\" font-size=\"14\" fill=\"white\">$(b.name)</text>")
      println(io, "<text x=\"$(r1(cx(b.index)))\" y=\"$(r1(y + 37))\" text-anchor=\"middle\" font-size=\"12\" fill=\"$(text)\">$(b.k == 0 ? "Slack" : string(b.type, ", k=", b.k))</text>")
    end
    base = H - 66
    println(io, "<text x=\"$(margin)\" y=\"$(base)\" font-size=\"12\" fill=\"#555550\">rows of J: 2k-1 (dP), 2k (dQ or dV); columns: k (Vr), $(length(structure.buses) - 1)+k (Vi)</text>")
    println(io, "<text x=\"$(margin)\" y=\"$(base + 20)\" font-size=\"12\" fill=\"#555550\">two circles: transformer; two lines: two systems; yellow frame: a bus of the example</text>")
    println(io, "<text x=\"$(margin)\" y=\"$(base + 40)\" font-size=\"12\" fill=\"#555550\">columns from left to right: distance from the slack bus in branches</text>")
    println(io, "</svg>")
  end
  return path
end

function print_block(title::String, b)
  println(title, "  (rows ", b.rows, ", columns ", b.cols, ")")
  @printf("    dP/dVr      = %12.6f    dP/dVi      = %12.6f\n", b.block[1, 1], b.block[1, 2])
  @printf("    d(Q|V)/dVr  = %12.6f    d(Q|V)/dVi  = %12.6f\n", b.block[2, 1], b.block[2, 2])
end

"""
    main(; casefile, eqBus, varBus, farBus)

Solves the shipped case and shows its structure first (console and
`network.svg`: bus types, the position k of every bus in the Jacobian,
lines, double systems, transformers). Then it builds the Jacobian at the
solution and reads the block of `eqBus` by `varBus` (two buses a line connects) and the one of
the opposite direction (four other values: the Jacobian is not
symmetric). The block of `eqBus` by `farBus`, a bus no branch connects it
to, is zero in every entry: the Jacobian has the sparsity of the
admittance matrix. Every block is checked against central differences of
the mismatch. Then the Newton iterations
from a flat start are repeated with the block recorded at every iterate
and printed as a table.
"""
function main(; casefile::String = joinpath(dirname(dirname(@__DIR__)), "data", "scf", "sp_case14.scf.json"), eqBus::String = "Roggau_110", varBus::String = "Ilmrode_110", farBus::String = "Hassel_20")
  print_example_banner("examples/others/exp_jacobian_block.jl", "the 2x2 block of the rectangular Jacobian between two buses, checked and drawn over the Newton iterations")

  net = Sparlectra.importSCF(casefile)
  ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg == 0 || error("power flow did not converge (erg = $(erg))")
  println("solved ", basename(casefile), " in ", ite, " iteration(s), ", length(net.nodeVec), " buses")
  println()

  p = jacobian_inputs(net)
  structure = network_structure(net, p.slack_idx)
  print_network(structure, [eqBus, varBus, farBus])
  diagram = write_network_diagram(joinpath(OUT_DIR, "network.svg"), structure, p.slack_idx, [eqBus, varBus, farBus])
  println("network diagram: ", relpath(diagram, dirname(dirname(@__DIR__))))
  println()
  J = jacobian_at(p, p.V)
  println("Jacobian at the solution: ", size(J, 1), " x ", size(J, 2), ", ", nnz(J), " stored entries")
  println()
  for (from, to) in ((eqBus, varBus), (varBus, eqBus), (eqBus, farBus))
    b = jacobianBlock(J, net, from, to; slack_idx = p.slack_idx)
    print_block("equations of $(from) by the voltage of $(to)", b)
    (to == farBus && iszero(b.block)) && println("    no branch connects the two buses: the block is structurally zero")
    fd = finite_difference_block(p, net, from, to)
    deviation = maximum(abs, fd .- b.block)
    @printf("    largest deviation from central differences: %.2e\n", deviation)
    deviation < 1e-6 || error("the block of $(from) by $(to) does not match the finite differences")
    println()
  end

  history = block_history(p, net, eqBus, varBus)
  println("the block of ", eqBus, " by ", varBus, " over the Newton iterations from a flat start")
  println("  it    max mismatch      dP/dVr      dP/dVi  d(Q|V)/dVr  d(Q|V)/dVi")
  for h in history
    @printf("  %2d    %12.3e  %10.5f  %10.5f  %10.5f  %10.5f\n", h.iteration, h.mismatch, h.block[1, 1], h.block[1, 2], h.block[2, 1], h.block[2, 2])
  end
  println()
  println("Reading the numbers:")
  println("  a block is nonzero only between buses a branch connects, and on the diagonal;")
  println("  the entries move most in the first iteration and settle with the voltages,")
  println("  which is what makes reusing a factorization over late iterations cheap.")
  return (block = jacobianBlock(J, net, eqBus, varBus; slack_idx = p.slack_idx), history = history, diagram = diagram)
end

run_example(main)
