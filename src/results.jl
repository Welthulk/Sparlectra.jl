# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 07.09.2023
# file results.jl
# Purpose: functions for formatting and printing results of power flow calculations
function format_version(version::VersionNumber)
  major = lpad(version.major, 2, '0')
  minor = lpad(version.minor, 1, '0')
  patch = lpad(version.patch, 2, '0')
  return "$major.$minor.$patch"
end

function formatBranchResults(net::Net)
  #! format: off
  formatted_results = @sprintf("\n===========================================================================================================================\n")

  formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Branch", "From", "To", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]")
  formatted_results *= @sprintf("===========================================================================================================================\n")
  #! format: on
  for br in net.branchVec
    from = br.fromBus
    to = br.toBus
    bName = br.comp.cName
    if br.status == 0 || (isnothing(br.fBranchFlow)) || (isnothing(br.tBranchFlow))
      pfromVal = qfromVal = ptoVal = qtoVal = pLossval = qLossval = 0.0
    else
      pfromVal = (br.fBranchFlow.pFlow === nothing) ? NaN : br.fBranchFlow.pFlow
      qfromVal = (br.fBranchFlow.qFlow === nothing) ? NaN : br.fBranchFlow.qFlow

      ptoVal = (br.tBranchFlow.pFlow === nothing) ? NaN : br.tBranchFlow.pFlow
      qtoVal = (br.tBranchFlow.qFlow === nothing) ? NaN : br.tBranchFlow.qFlow

      pLossval = (br.pLosses === nothing) ? NaN : br.pLosses
      qLossval = (br.qLosses === nothing) ? NaN : br.qLosses
      ratedS = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA

      check = false
      if ratedS > 0.0
        if max(abs(pfromVal), abs(ptoVal)) > ratedS
          check = true
        end
      end

      if check
        bName *= " !"
      end
    end
    #! format: off
    formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f |  %-10.3f|\n", bName, from, to, pfromVal, qfromVal, ptoVal, qtoVal, pLossval, qLossval)
    #! format: on
  end
  formatted_results *= @sprintf("---------------------------------------------------------------------------------------------------------------------------\n")
  (∑pv, ∑qv) = getTotalLosses(net = net)
  total_losses = @sprintf("total network power balance (Σ S_branch): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

  return formatted_results, total_losses
end

function printACPFlowResults(net::Net, ct::Float64, ite::Int, tol::Float64, toFile::Bool = false, path::String = ""; converged::Bool = true, solver::Symbol = :NR)
  if toFile
    filename = strip("result_$(net.name).txt")
    io = open(joinpath(path, filename), "w")
    @info "Results are written to $(joinpath(path, filename))"
  else
    io = Base.stdout
  end

  vers = Sparlectra.SparlectraVersion
  current_date = Dates.format(Dates.now(), "d-u-yy H:M:S")

  formatted_version = format_version(vers)
  flowResults, totalLosses = formatBranchResults(net)

  @printf(io, "================================================================================\n")
  @printf(io, "| SPARLECTRA Version %-10s - AC Power Flow Results                        |\n", formatted_version)
  @printf(io, "================================================================================\n")

  busses = length(net.nodeVec)
  branches = length(net.branchVec)
  lines = length(net.linesAC)
  trafos = length(net.trafos)
  gens = 0
  loads = 0
  shunts = 0
  auxb = 0

  nodes = sort(net.nodeVec, by = x -> x.busIdx)

  npv = 0
  npq = 0
  niso = 0
  for n in nodes
    npv += n._nodeType == Sparlectra.PV ? 1 : 0
    npq += isPQNode(n) ? 1 : 0
    niso += isIsolated(n) ? 1 : 0
    if occursin("_Aux_", n.comp.cName)
      auxb += 1
    end
  end
  for ps in net.prosumpsVec
    loads += ps.proSumptionType == Sparlectra.Consumption ? 1 : 0
    gens += ps.proSumptionType == Sparlectra.Injection ? 1 : 0
  end
  shunts = length(net.shuntVec)

  @printf(io, "Date           :%20s\n", current_date)
  @printf(io, "Iterations     :%10d\n", ite)
  @printf(io, "Flatstart      :%10s\n", net.flatstart ? "Yes" : "No")
  @printf(io, "Tolerance      : %.1e\n", tol)
  @printf(io, "Solver         :%15s\n", string(solver))
  if converged
    @printf(io, "Converged in   :%10f seconds\n", ct)
  else
    @printf(io, "Status         :%10s\n", "Not Converged")
  end
  @printf(io, "Case           :%15s\n", net.name)
  @printf(io, "BaseMVA        :%10d\n", net.baseMVA)
  if auxb > 0 && niso > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Iso: %d Slack: %d\n", busses, npv, npq, auxb, niso, 1)
  elseif auxb > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Slack: %d\n", busses, npv, npq, auxb, 1)
  else
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d Slack: %d)\n", busses, npv, npq, 1)
  end
  @printf(io, "Branches       :%10d\n", branches)
  @printf(io, "Lines          :%10d\n", lines)
  @printf(io, "Trafos         :%10d\n", trafos)
  @printf(io, "Generators     :%10d\n", gens)
  @printf(io, "Loads          :%10d\n", loads)
  @printf(io, "Shunts         :%10d\n", shunts)

  num_q_limit = length(net.qLimitEvents)
  @printf(io, "PV→PQ (Q-Limit):%10d\n", num_q_limit)

  println(io, "\n", totalLosses)

  @printf(io, "===============================================================================================================================================================================\n")
  @printf(io, "| %-5s | %-20s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Nr", "Bus", "Vn [kV]", "V [kV]", "V [pu]", "phi [deg]", "Pg [MW]", "Qg [MVar]", "Pl [MW]", "Ql [MVar]", "Ps [MW]", "Qs [MVar]", "Type")
  @printf(io, "===============================================================================================================================================================================\n")

  pGS = qGS = pLS = qLS = ""
  tpGS = tqGS = tpLS = tqLS = 0.0
  pShunt_str = qShunt_str = ""
  tpShunt = tqShunt = 0.0
  for n in nodes
    if !isnothing(n._pƩGen)
      abs(n._pƩGen) > 1e-6
      pGS = @sprintf("%10.3f", n._pƩGen)
      tpGS += n._pƩGen
    else
      pGS = ""
    end
    if !isnothing(n._qƩGen) && abs(n._qƩGen) > 1e-6
      qGS = @sprintf("%10.3f", n._qƩGen)
      tqGS += n._qƩGen
    else
      qGS = ""
    end

    if !isnothing(n._pƩLoad) && abs(n._pƩLoad) > 1e-6
      pLS = @sprintf("%10.3f", n._pƩLoad)
      tpLS += n._pƩLoad
    else
      pLS = ""
    end
    if !isnothing(n._qƩLoad) && abs(n._qƩLoad) > 1e-6
      qLS = @sprintf("%10.3f", n._qƩLoad)
      tqLS += n._qƩLoad
    else
      qLS = ""
    end
    if !isnothing(n._pShunt) && abs(n._pShunt) > 1e-6
      pShunt_str = @sprintf("%10.3f", n._pShunt)
      tpShunt += n._pShunt
    else
      pShunt_str = ""
    end
    if !isnothing(n._qShunt) && abs(n._qShunt) > 1e-6
      qShunt_str = @sprintf("%10.3f", n._qShunt)
      tqShunt += n._qShunt
    else
      qShunt_str = ""
    end
    typeStr = toString(n._nodeType)

    # Mark PV→PQ buses (hit Q-limit) with a star in the Type column
    if haskey(net.qLimitEvents, n.busIdx)
      typeStr *= "*"
    end

    v = n.comp.cVN * n._vm_pu
    nodeName = n.comp.cName
    if !isnothing(n._vmin_pu) && !isnothing(n._vmax_pu)
      if !isIsolated(n) && (n._vm_pu < n._vmin_pu || n._vm_pu > n._vmax_pu)
        nodeName *= " !"
      end
    end

    @printf(io, "| %-5d | %-20s | %-10.1f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", n.busIdx, nodeName, n.comp.cVN, v, n._vm_pu, n._va_deg, pGS, qGS, pLS, qLS, pShunt_str, qShunt_str, typeStr)
  end

  @printf(io, "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  println(io, flowResults)

  if toFile
    close(io)
    #println("Results have been written to $(joinpath(path, filename))")
  end
end

function formatProsumerResults(net::Net)
  buf = IOBuffer()

  # Rebuild mapping: bus index -> bus name
  busNameByIdx = Dict{Int,String}()
  for (name, idx) in net.busDict
    busNameByIdx[idx] = name
  end

  # Collect indices per bus, separated into generators and loads
  gens_by_bus  = Dict{Int,Vector{Int}}()
  loads_by_bus = Dict{Int,Vector{Int}}()

  for (idx, ps) in enumerate(net.prosumpsVec)
    bus = getPosumerBusIndex(ps)
    if isGenerator(ps)
      push!(get!(gens_by_bus, bus, Int[]), idx)
    else
      push!(get!(loads_by_bus, bus, Int[]), idx)
    end
  end

  # Helper: compute generator status from Q and minQ/maxQ
  status_from_Q = function (ps, qres)
    # default status
    status = "ok"
    isnothing(qres) && return status

    # tolerance for limit detection
    tol = 1e-6

    if !isnothing(ps.maxQ) && qres >= ps.maxQ - tol
      status = "Q-max limit"
    elseif !isnothing(ps.minQ) && qres <= ps.minQ + tol
      status = "Q-min limit"
    end

    return status
  end

  # =========================
  # Generator section
  # =========================
  println(buf, "\nGenerator results:")
  println(buf, "────────────────────────────────────────────────────────")
  @printf(buf, "%-8s %4s %12s %12s   %-14s\n", "Bus", "Gen#", "P [MW]", "Q [MVar]", "Status")
  println(buf, "────────────────────────────────────────────────────────")

  for bus in sort(collect(keys(gens_by_bus)))
    gens_idx = gens_by_bus[bus]

    # optional: sort generators at a bus by component name
    sort!(gens_idx, by = i -> net.prosumpsVec[i].comp.cName)

    busName = get(busNameByIdx, bus, "Bus_$bus")

    for (k, idx) in enumerate(gens_idx)
      ps = net.prosumpsVec[idx]

      # Prefer results (pRes/qRes); fall back to spec values
      p = ps.pRes === nothing ? ps.pVal : ps.pRes
      q = ps.qRes === nothing ? ps.qVal : ps.qRes

      status = status_from_Q(ps, q)

      @printf(buf, "%-8s %4d %12.3f %12.3f   %-14s\n", busName, k, p === nothing ? 0.0 : p, q === nothing ? 0.0 : q, status)
    end
  end

  println(buf, "────────────────────────────────────────────────────────")

  # =========================
  # Load section
  # =========================
  println(buf, "\nLoad results:")
  println(buf, "────────────────────────────────────────")
  @printf(buf, "%-8s %5s %12s %12s\n", "Bus", "Load#", "P [MW]", "Q [MVar]")
  println(buf, "────────────────────────────────────────")

  if !isnothing(loads_by_bus)
    for bus in sort(collect(keys(loads_by_bus)))
      try
        loads_idx = loads_by_bus[bus]

        # optional: sort loads at a bus by component name
        sort!(loads_idx, by = i -> net.prosumpsVec[i].comp.cName)

        busName = get(busNameByIdx, bus, "Bus_$bus")

        for (k, idx) in enumerate(loads_idx)
          ps = net.prosumpsVec[idx]

          # Prefer results (pRes/qRes); fall back to spec values
          p = ps.pRes === nothing ? ps.pVal : ps.pRes
          q = ps.qRes === nothing ? ps.qVal : ps.qRes

          @printf(buf, "%-8s %5d %12.3f %12.3f\n", busName, k, p === nothing ? 0.0 : p, q === nothing ? 0.0 : q)
        end
      catch e
        @warn "Error formatting load results for bus index $bus: $e"
      end
    end
  end

  println(buf, "────────────────────────────────────────────────────────")

  return String(take!(buf))
end

function printProsumerResults(net::Net)
  prosText = formatProsumerResults(net)
  println(prosText)
end
