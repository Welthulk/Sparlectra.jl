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

# Advanced Chaos Analysis with Visualization and Detailed Investigation
# Goal: Visualization of Basin of Attraction and critical regions

using Sparlectra
using Logging
using Printf
global_logger(ConsoleLogger(stderr, Logging.Info))


function converges_for(r_factor, l_factor; v0=1.0, a0=0.0, max_iter=50, tol=1e-6,
                       method=:rectangular, opt_sparse=true, verbose=0)
  res = test_convergence(
    voltage_B2 = v0,
    angle_B2   = a0,
    angle_B3   = 0.0,
    resistance_factor = r_factor,
    load_factor       = l_factor,
    impedance_factor  = r_factor,
    max_iter = max_iter,
    tolerance = tol,
    method = method,
    opt_sparse = opt_sparse,
    verbose = verbose,
  )
  return res.converged, res.iterations
end

# --- Drop-in: adaptive critical-load search (expand upper bound + bisection) ---
# Put this block near your existing helper functions.

"""
    find_load_critical_adaptive(r_factor;
        l_lo=1.0, l_hi=6.0,
        grow=2.0, l_max=50.0,
        itmax=30, eps=1e-3,
        kwargs...)

Find approximate critical load factor L_crit for a given r_factor.
Strategy:
1) Ensure l_lo converges (otherwise return immediately: "already diverged at lower bound")
2) Expand l_hi geometrically until divergence is observed or l_max reached
3) Bisection on [l_ok, l_fail] until eps

Returns a NamedTuple: (lcrit, bracket=(a,b), note)
- lcrit is the last known convergent value (≈ lower end of boundary)
- bracket is (l_ok, l_fail) if divergence found, else (l_lo, l_hi)
"""
function find_load_critical_adaptive(r_factor;
    l_lo::Float64 = 1.0,
    l_hi::Float64 = 6.0,
    grow::Float64 = 2.0,
    l_max::Float64 = 50.0,
    itmax::Int = 30,
    eps::Float64 = 1e-3,
    kwargs...,
)
    ok_lo, _ = converges_for(r_factor, l_lo; kwargs...)
    if !ok_lo
        return (lcrit = l_lo, bracket = (l_lo, l_hi), note = "already diverged at lower bound")
    end

    # Expand upper bound until it fails (or we hit l_max)
    a = l_lo
    b = l_hi
    ok_b, _ = converges_for(r_factor, b; kwargs...)

    while ok_b && b < l_max
        a = b
        b = min(l_max, b * grow)
        ok_b, _ = converges_for(r_factor, b; kwargs...)
        # If b == l_max and still ok, we're done
        if b >= l_max && ok_b
            return (lcrit = b, bracket = (l_lo, b), note = "still converged at l_max")
        end
    end

    # If it fails already at initial l_hi, we already have a bracket [a,b] with a ok, b fail
    # If it only failed after expansion, also bracketed.
    # Now bisection on [a,b] (a converged, b diverged)
    for _ = 1:itmax
        m = 0.5 * (a + b)
        ok_m, _ = converges_for(r_factor, m; kwargs...)
        if ok_m
            a = m
        else
            b = m
        end
        (b - a) < eps && break
    end

    return (lcrit = a, bracket = (a, b), note = "adaptive+bisection")
end

function find_load_critical(r_factor; l_lo=1.0, l_hi=6.0, itmax=25, eps=1e-3, kwargs...)
  ok_lo, _ = converges_for(r_factor, l_lo; kwargs...)
  ok_hi, _ = converges_for(r_factor, l_hi; kwargs...)

  # Falls hi noch konvergiert: Grenze liegt höher als l_hi
  if ok_hi
    return (lcrit = l_hi, bracket = (l_lo, l_hi), note = "still converged at upper bound")
  end
  # Falls lo schon nicht konvergiert: Grenze liegt unter l_lo
  if !ok_lo
    return (lcrit = l_lo, bracket = (l_lo, l_hi), note = "already diverged at lower bound")
  end

  a = l_lo
  b = l_hi
  for _ = 1:itmax
    m = 0.5*(a+b)
    ok, _ = converges_for(r_factor, m; kwargs...)
    if ok
      a = m
    else
      b = m
    end
    (b-a) < eps && break
  end
  return (lcrit = a, bracket = (a, b), note = "bisection")
end

# Helper functions for safe network operations
function create_base_network()
  net = Net(name = "BaseNetwork", baseMVA = 100.0, vmin_pu = 0.1, vmax_pu = 2.0)
  base_r, base_x = 0.02, 0.08

  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", busType = "PV", vn_kV = 110.0)

  
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 15.0, r =  base_r, x =  base_x, b = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 20.0, r =  (base_r + 0.01), x =  (base_x + 0.02), b = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r =  (base_r + 0.02), x =  (base_x + 0.04), b = 0.0, ratedS = 100.0)

  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 30.0)
  addProsumer!(net=net, busName="B3", type="SYNCHRONOUSMACHINE",  p=50.0, vm_pu=1.05, qMin=-300.0, qMax=300.0)  
  addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION",  vm_pu=1.0, va_deg=0.0, referencePri="B1")
  return net
end

function create_variable_network(; resistance_factor = 1.0, load_factor = 1.0, impedance_factor = 1.0)
  net = Net(name = "VarNet_R$(resistance_factor)_L$(load_factor)_X$(impedance_factor)", baseMVA = 100.0, vmin_pu = 0.1, vmax_pu = 2.0)

  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", busType = "PV", vn_kV = 110.0)

  base_r, base_x = 0.02, 0.08
  r_mult = resistance_factor
  x_mult = impedance_factor

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 15.0, r =  base_r * r_mult, x =  base_x * x_mult, b = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 20.0, r =  (base_r + 0.01) * r_mult, x =  (base_x + 0.02) * x_mult, b = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r =  (base_r + 0.02) * r_mult, x =  (base_x + 0.04) * x_mult, b = 0.0, ratedS = 100.0)

  base_p, base_q = 40.0, 30.0
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = base_p * load_factor, q = base_q * load_factor)
  addProsumer!(net=net, busName="B3", type="SYNCHRONOUSMACHINE",  p=50.0, vm_pu=1.05, qMin=-300.0, qMax=300.0)  
  addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION",  vm_pu=1.0, va_deg=0.0, referencePri="B1")

  return net
end

function safe_set_initial_values!(net; voltage_B2 = 1.0, angle_B2 = 0.0, angle_B3 = 0.0)
  success = true
  try
    if isdefined(Main, :setNodeVoltage!)
      setNodeVoltage!(net = net, busName = "B2", vm_pu = voltage_B2, va_deg = angle_B2)
    end
    if isdefined(Main, :setNodeAngle!)
      setNodeAngle!(net = net, busName = "B3", va_deg = angle_B3)
    end
  catch e
    success = false
  end
  return success
end

function test_convergence(;
    voltage_B2 = 1.0,
    angle_B2   = 0.0,
    angle_B3   = 0.0,
    resistance_factor = 1.0,
    load_factor       = 1.0,
    impedance_factor  = 1.0,
    max_iter::Int = 50,
    tolerance::Float64 = 1e-6,
    method::Symbol = :rectangular,
    opt_sparse::Bool = true,
    verbose::Int = 0,
)
  try
    net = create_variable_network(
      resistance_factor = resistance_factor,
      load_factor       = load_factor,
      impedance_factor  = impedance_factor,
    )

    ok, msg = validate!(net = net)
    ok || return (converged = false, iterations = 0)

    safe_set_initial_values!(net,
      voltage_B2 = float(voltage_B2),
      angle_B2   = float(angle_B2),
      angle_B3   = float(angle_B3),
    )

    iterations, result = runpf!(
      net, max_iter, tolerance, verbose;
      method = method,
      opt_sparse = opt_sparse,
    )

    return (converged = (result == 0), iterations = iterations)
  catch
    return (converged = false, iterations = 0)
  end
end

# 1. Basin of Attraction Analysis for Starting Values
function basin_of_attraction_analysis()
  println("🔬 BASIN OF ATTRACTION ANALYSIS - STARTING VALUES")
  println(repeat("=", 60))

  # Define resolution for analysis
  voltage_range = 0.4:0.1:1.6
  angle_range = -20:5:25

  println("Analyzing $(length(voltage_range)) × $(length(angle_range)) = $(length(voltage_range) * length(angle_range)) starting value combinations")
  println("Voltage range: $(first(voltage_range)) - $(last(voltage_range)) pu")
  println("Angle range: $(first(angle_range))° - $(last(angle_range))°")
  println()

  convergence_map = zeros(Int, length(voltage_range), length(angle_range))
  iteration_map = zeros(Int, length(voltage_range), length(angle_range))

  total_tests = length(voltage_range) * length(angle_range)
  current_test = 0

  println("Progress: [" * repeat(" ", 50) * "]")
  print("Progress: [")

  for (i, v_start) in enumerate(voltage_range)
    for (j, a_start) in enumerate(angle_range)
      current_test += 1

      # Progress Bar
      if current_test % max(1, div(total_tests, 50)) == 0
        print("▓")
      end

      result = test_convergence(voltage_B2 = v_start, angle_B2 = a_start, angle_B3 = 0.0, resistance_factor = 1.0, load_factor = 1.0, impedance_factor = 1.0)

      convergence_map[i, j] = result.converged ? 1 : 0
      iteration_map[i, j] = result.iterations
    end
  end

  println("]")
  println()

  # Analysis of results
  total_points = length(convergence_map)
  converged_points = sum(convergence_map)
  chaos_rate = round(100 * (total_points - converged_points) / total_points, digits = 1)

  println("📊 BASIN OF ATTRACTION RESULTS:")
  println("Tested points: $total_points")
  println("Converged points: $converged_points")
  println("Chaos rate: $chaos_rate%")
  println()

  # ASCII Visualization
  println("🎨 BASIN OF ATTRACTION MAP:")
  println("Legend: ✓ = Converged, ✗ = Diverged")
  println()
  print("V\\θ ")
  for angle in angle_range
    @printf("%4d°", angle)
  end
  println()

  for (i, voltage) in enumerate(voltage_range)
    @printf("%.1f", voltage)
    for j = 1:length(angle_range)
      if convergence_map[i, j] == 1
        print("  ✓ ")
      else
        print("  ✗ ")
      end
    end
    println()
  end

  # Detailed iteration map
  println()
  println("🔢 ITERATION MAP:")
  println("Shows number of iterations until convergence (0 = Diverged)")
  println()
  print("V\\θ ")
  for angle in angle_range
    @printf("%4d°", angle)
  end
  println()

  for (i, voltage) in enumerate(voltage_range)
    @printf("%.1f", voltage)
    for j = 1:length(angle_range)
      @printf("%4d ", iteration_map[i, j])
    end
    println()
  end

  return convergence_map, iteration_map, collect(voltage_range), collect(angle_range)
end

# 2. Parameter Sensitivity Analysis
function parameter_sensitivity_analysis()
  println("\n⚙️ PARAMETER SENSITIVITY ANALYSIS (critical load boundary, adaptive)")
  println(repeat("=", 72))

  resistance_range = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]

  println("R-factor | L_crit (approx) | bracket            | note")
  println(repeat("-", 72))

  for r in resistance_range
    res = find_load_critical_adaptive(r;
      l_lo=1.0, l_hi=6.0, grow=2.0, l_max=50.0, eps=1e-3, itmax=30,
      method=:rectangular, opt_sparse=true, verbose=0,
      max_iter=50, tol=1e-6,
    )
    @printf("%7.1fx | %14.3f | [%-7.3f, %-7.3f] | %s\n",
      r, res.lcrit, res.bracket[1], res.bracket[2], res.note)
  end

  # Keep return shape compatible with your existing call site (expects 4 things).
  return nothing, nothing, resistance_range, nothing
end

# 3. Critical Boundary Analysis
function critical_boundary_analysis()
  println("\n🎯 CRITICAL BOUNDARY ANALYSIS")
  println(repeat("=", 60))

  println("Searching for exact transitions between convergence and divergence...")

  # Find critical voltage boundary
  println("\n🔍 Critical voltage boundary (at θ=0°):")
  critical_voltages = []

  for voltage = 0.3:0.05:0.6
    result = test_convergence(voltage_B2 = voltage, angle_B2 = 0.0, angle_B3 = 0.0)
    status = result.converged ? "✓" : "✗"
    iterations = result.iterations

    @printf("V = %.2f pu: %s (%d Iter.)\n", voltage, status, iterations)

    if !result.converged
      push!(critical_voltages, voltage)
    end
  end

  # Find critical angle boundary
  println("\n🔍 Critical angle boundary (at V=1.0 pu):")
  critical_angles = []

  for angle = 20:5:45
    result = test_convergence(voltage_B2 = 1.0, angle_B2 = angle, angle_B3 = 0.0)
    status = result.converged ? "✓" : "✗"
    iterations = result.iterations

    @printf("θ = %d°: %s (%d Iter.)\n", angle, status, iterations)

    if !result.converged
      push!(critical_angles, angle)
    end
  end

  # Find critical parameter combinations
  println("\n🔍 Critical parameter combinations:")
  critical_params = []

  for r_factor in [3.0, 5.0, 8.0, 10.0]
    for l_factor in [2.0, 3.0, 4.0]
      result = test_convergence(voltage_B2 = 1.0, angle_B2 = 0.0, angle_B3 = 0.0, resistance_factor = r_factor, load_factor = l_factor, impedance_factor = r_factor)
      status = result.converged ? "✓" : "✗"
      iterations = result.iterations

      @printf("R=%.1fx, L=%.1fx: %s (%d Iter.)\n", r_factor, l_factor, status, iterations)

      if !result.converged
        push!(critical_params, (r_factor, l_factor))
      end
    end
  end

  println("\n📋 CRITICAL REGIONS SUMMARY:")
  if !isempty(critical_voltages)
    println("Critical voltages: $(minimum(critical_voltages)) - $(maximum(critical_voltages)) pu")
  end
  if !isempty(critical_angles)
    println("Critical angles: $(minimum(critical_angles))° - $(maximum(critical_angles))°")
  end
  if !isempty(critical_params)
    println("Critical parameter combinations found: $(length(critical_params))")
  end

  return critical_voltages, critical_angles, critical_params
end

# 4. Combined 3D Analysis
function combined_3d_analysis()
  println("\n🌪️ COMBINED 3D CHAOS ANALYSIS")
  println(repeat("=", 60))

  println("Analyzing interaction between starting values AND parameters...")

  # Reduced resolution for 3D analysis
  voltage_range = [0.6, 0.8, 1.0, 1.2, 1.4]
  angle_range = [-15, -5, 0, 5, 15]
  param_range = [(1.0, 1.0), (2.0, 1.5), (3.0, 2.0), (5.0, 2.5), (8.0, 3.0)]

  chaos_scenarios = []

  println("Testing $(length(voltage_range))×$(length(angle_range))×$(length(param_range)) = $(length(voltage_range)*length(angle_range)*length(param_range)) combinations")
  println()

  scenario_count = 0
  chaos_count = 0

  for v_start in voltage_range
    for a_start in angle_range
      for (r_factor, l_factor) in param_range
        scenario_count += 1

        result = test_convergence(voltage_B2 = v_start, angle_B2 = a_start, angle_B3 = 0.0, resistance_factor = r_factor, load_factor = l_factor, impedance_factor = r_factor)

        if !result.converged
          chaos_count += 1
          push!(chaos_scenarios, (v_start, a_start, r_factor, l_factor, result.iterations))
        end
      end
    end
  end

  combined_chaos_rate = round(100 * chaos_count / scenario_count, digits = 1)

  println("📊 COMBINED 3D ANALYSIS RESULTS:")
  println("Tested scenarios: $scenario_count")
  println("Chaos scenarios: $chaos_count")
  println("3D chaos rate: $combined_chaos_rate%")
  println()

  if !isempty(chaos_scenarios)
    println("🎯 TOP-10 CHAOS HOTSPOTS:")
    # Sort by "extremeness" (low voltage + high parameters)
    sorted_chaos = sort(chaos_scenarios, by = x -> x[1] + x[3] + x[4], rev = false)

    for i = 1:min(10, length(sorted_chaos))
      v, a, r, l, iter = sorted_chaos[i]
      @printf("%2d. V=%.1f pu, θ=%3d°, R=%.1fx, L=%.1fx (Iter: %d)\n", i, v, Int(a), r, l, iter)
    end
  end

  return chaos_scenarios, combined_chaos_rate
end

# 5. Statistical Chaos Report
function generate_chaos_report()
  println("\n📊 STATISTICAL CHAOS REPORT")
  println(repeat("=", 60))

  # Collect data for different categories
  categories = [
    ("StartVal_Mild", [(1.0, 0.0, 1.0, 1.0), (0.9, -5.0, 1.0, 1.0), (1.1, 5.0, 1.0, 1.0)]),
    ("StartVal_Extreme", [(0.5, -20.0, 1.0, 1.0), (1.8, 25.0, 1.0, 1.0), (0.3, 35.0, 1.0, 1.0)]),
    ("Param_Mild", [(1.0, 0.0, 2.0, 1.5), (1.0, 0.0, 1.5, 2.0), (1.0, 0.0, 2.5, 1.0)]),
    ("Param_Extreme", [(1.0, 0.0, 8.0, 3.0), (1.0, 0.0, 10.0, 4.0), (1.0, 0.0, 15.0, 5.0)]),
    ("Combined_Mild", [(0.8, -10.0, 2.0, 2.0), (1.2, 10.0, 2.0, 2.0), (0.9, -5.0, 3.0, 1.5)]),
    ("Combined_Extreme", [(0.4, -25.0, 8.0, 4.0), (1.6, 30.0, 10.0, 5.0), (0.2, 40.0, 15.0, 6.0)]),
  ]

  println("Category             | Chaos Rate | Avg. Iter. | Stability")
  println(repeat("-", 60))

  for (category_name, scenarios) in categories
    chaos_count = 0
    total_iterations = 0

    for (v, a, r, l) in scenarios
      result = test_convergence(voltage_B2 = v, angle_B2 = a, angle_B3 = 0.0, resistance_factor = r, load_factor = l, impedance_factor = r)

      if !result.converged
        chaos_count += 1
      end
      total_iterations += result.iterations
    end

    chaos_rate = round(100 * chaos_count / length(scenarios), digits = 1)
    avg_iterations = round(total_iterations / length(scenarios), digits = 1)
    stability = chaos_rate < 30 ? "Stable" : chaos_rate < 70 ? "Moderate" : "Unstable"

    @printf("%-20s | %7.1f%% | %8.1f | %-8s\n", category_name, chaos_rate, avg_iterations, stability)
  end

  println()
  println("🎯 CHAOS CHARACTERISTICS:")
  println("• Starting values: Local instabilities")
  println("• Parameters: Global system weakness")
  println("• Combination: Enhanced nonlinearity")
  println("• Iterations: Convergence speed as indicator")
end

# Main function: Comprehensive Chaos Visualization
function comprehensive_chaos_visualization()
  println("🎨 COMPREHENSIVE CHAOS VISUALIZATION AND ANALYSIS")
  println(repeat("=", 70))
  println("Performing detailed investigation of all chaos aspects")
  println()

  # 1. Basin of Attraction
  basin_conv, basin_iter, v_range, a_range = basin_of_attraction_analysis()

  # 2. Parameter Sensitivity  
  param_conv, param_iter, r_range, l_range = parameter_sensitivity_analysis()

  # 3. Critical Boundaries
  crit_v, crit_a, crit_p = critical_boundary_analysis()

  # 4. 3D Combination
  chaos_hotspots, combined_rate = combined_3d_analysis()

  # 5. Statistical Report
  generate_chaos_report()

  println("\n" * repeat("=", 70))
  println("🏆 CHAOS VISUALIZATION COMPLETED")
  println("All critical regions and chaos patterns identified!")

  return (basin_data = (basin_conv, basin_iter, v_range, a_range), param_data = (param_conv, param_iter, r_range, l_range), critical_data = (crit_v, crit_a, crit_p), hotspots = chaos_hotspots, combined_chaos_rate = combined_rate)
end

function parameter_sensitivity_plot_ascii(; method=:rectangular, opt_sparse=true, verbose=0,
                                          max_iter=50, tol=1e-6)
  resistance_range = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]

  # collect results
  R = Float64[]
  L = Float64[]
  br_lo = Float64[]
  br_hi = Float64[]
  note = String[]

  for r in resistance_range
    res = find_load_critical(r;
      l_lo=1.0, l_hi=2.0,   # start low; adaptive will grow
      eps=1e-3, itmax=25,
      method=method, opt_sparse=opt_sparse, verbose=verbose,
      max_iter=max_iter, tol=tol,
    )
    push!(R, r)
    push!(L, res.lcrit)
    push!(br_lo, res.bracket[1])
    push!(br_hi, res.bracket[2])
    push!(note, res.note)
  end

  println("\n📈 Lcrit vs R (ASCII)")
  println("R-factor | L_crit | bracket                | note")
  println(repeat("-", 72))
  for i in eachindex(R)
    @printf("%7.1fx | %6.3f | [%-6.3f, %-6.3f] | %s\n",
            R[i], L[i], br_lo[i], br_hi[i], note[i])
  end

  # ASCII bar plot
  # scale bars to max Lcrit
  Lmax = maximum(L)
  width = 50
  println("\nASCII plot (bar length ~ Lcrit):")
  for i in eachindex(R)
    n = Int(round(width * (L[i] / Lmax)))
    bar = repeat("█", n)
    @printf("R=%4.1fx | %-50s  %.3f\n", R[i], bar, L[i])
  end

  return (R=R, L=L, br_lo=br_lo, br_hi=br_hi, note=note)
end

# === MAIN EXECUTION ===
println("🚀 Starting comprehensive chaos visualization and analysis...")
println()

# Perform complete analysis
results = comprehensive_chaos_visualization()
sens = parameter_sensitivity_plot_ascii()
println("\n🎉 CHAOS ANALYSIS COMPLETE!")
println("All Basin of Attraction, critical regions and chaos hotspots have been identified!")
println("This data can be used for scientific visualizations.")
