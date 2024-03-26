# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 07.09.2023
# include-file results.jl
using Dates

function format_version(version::VersionNumber)
  major = lpad(version.major, 2, '0')
  minor = lpad(version.minor, 2, '0')
  patch = lpad(version.patch, 2, '0')
  return "$major.$minor.$patch"
end

function formatBranchResults(net::ResDataTypes.Net)

  formatted_results = @sprintf("\n===========================================================================================================================\n")

  formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Branch", "From", "To", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]")
  formatted_results *= @sprintf("===========================================================================================================================\n")

  ∑pv = ∑qv = 0.0
  ∑pfrom=∑qfrom=∑pto=∑qto=0.0
  
  for br in net.branchVec
    if br.status == 0
      continue
    end
    from = br.fromBus
    to = br.toBus
    bName = br.comp.cName
    pfromVal = (br.fBranchFlow.pFlow === nothing) ? NaN : br.fBranchFlow.pFlow
    qfromVal = (br.fBranchFlow.qFlow === nothing) ? NaN : br.fBranchFlow.qFlow

    ptoVal = (br.tBranchFlow.pFlow === nothing) ? NaN : br.tBranchFlow.pFlow
    qtoVal = (br.tBranchFlow.qFlow === nothing) ? NaN : br.tBranchFlow.qFlow

    VmFrom = (br.fBranchFlow.vm_pu === nothing) ? NaN : br.fBranchFlow.vm_pu
    VaFrom = (br.fBranchFlow.va_deg === nothing) ? NaN : br.fBranchFlow.va_deg

    VmTo = (br.tBranchFlow.vm_pu === nothing) ? NaN : br.tBranchFlow.vm_pu
    VaTo = (br.tBranchFlow.va_deg === nothing) ? NaN : br.tBranchFlow.va_deg
    ratedS = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
    ∑pfrom += pfromVal
    ∑qfrom += qfromVal
    ∑pto += ptoVal
    ∑qto += qtoVal

    v1 = VmFrom * exp(im * deg2rad(VaFrom))
    v2 = VmTo * exp(im * deg2rad(VaTo))
    # Losses: abs( Vf / tau - Vt ) ^ 2 / (Rs - j Xs)
    ratio = (br.ratio != 0.0) ? br.ratio : 1.0
    angle = (br.ratio != 0.0) ? br.angle : 0.0
    tap = calcComplexRatio(ratio, angle)
    

    Vdiff = v1 / tap - v2

    y = inv((br.r_pu + br.x_pu * im)) #+ 0.5*(br.g_pu + br.b_pu*im)

    Sdiff = Vdiff * conj(Vdiff * y) * net.baseMVA
    pv = real(Sdiff)
    qv = imag(Sdiff)
    ∑pv += abs(pv)
    ∑qv += abs(qv)
    check = false
    if ratedS > 0.0
      if max(abs(pfromVal), abs(ptoVal)) > ratedS
        check = true
      end
    end
    
    if check
      bName *= " !"
    end  
    formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f |  %-10.3f|\n", bName, from, to, pfromVal, qfromVal, ptoVal, qtoVal, abs(pv), abs(qv))
    
  end
  formatted_results *= @sprintf("---------------------------------------------------------------------------------------------------------------------------\n")  
  #formatted_results *= @sprintf("| %-10.3f | %-10.3f | %-10.3f | %-10.3f |\n", ∑pfrom, ∑qfrom, ∑pto, ∑qto)

  total_losses = @sprintf("total losses (I^2*Z): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

  return formatted_results, total_losses
end

function printACPFlowResults(net::ResDataTypes.Net, ct::Float64, ite::Int, tol::Float64, toFile::Bool = false, path::String = "")
  if toFile
    filename = strip("result_$(net.name).txt")
    io = open(joinpath(path, filename), "w")
    @info "Results are written to $(joinpath(path, filename))"
  else
    io = Base.stdout
  end

  vers = ResDataTypes.SparlectraVersion
  current_date = Dates.format(Dates.now(), "d-u-yy H:M:S")

  formatted_version = format_version(vers)
  flowResults, totalLosses = formatBranchResults(net)

  @printf(io,"================================================================================\n")
  @printf(io,"| SPARLECTRA Version %-10s - AC Power Flow Results                        |\n", formatted_version)
  @printf(io,"================================================================================\n")


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
  for n in nodes
    npv += n._nodeType == ResDataTypes.PV ? 1 : 0
    npq += n._nodeType == ResDataTypes.PQ ? 1 : 0
    if occursin("_Aux_", n.comp.cName)
      auxb += 1
    end
  end
  for ps in net.prosumpsVec
    loads += ps.proSumptionType == ResDataTypes.Consumption ? 1 : 0
    gens += ps.proSumptionType == ResDataTypes.Injection ? 1 : 0
  end
  shunts = length(net.shuntVec)

  @printf(io, "Date         :%20s\n", current_date)
  @printf(io, "Iterations   :%10d\n", ite)
  @printf(io, "Tolerance    : %.1e\n", tol) 
  @printf(io, "Converged in :%10f seconds\n", ct)
  @printf(io, "Case         :%15s\n", net.name)
  @printf(io, "BaseMVA      :%10d\n", net.baseMVA)
  if auxb > 0
    @printf(io, "Nodes        :%10d (PV: %d PQ: %d (Aux: %d) Slack: %d\n", busses, npv, npq, auxb, 1)
  else
    @printf(io, "Nodes        :%10d (PV: %d PQ: %d Slack: %d)\n", busses, npv, npq, 1)
  end
  @printf(io, "Branches     :%10d\n", branches)
  @printf(io, "Lines        :%10d\n", lines)
  @printf(io, "Trafos       :%10d\n", trafos)
  @printf(io, "Generators   :%10d\n", gens)
  @printf(io, "Loads        :%10d\n", loads)
  @printf(io, "Shunts       :%10d\n", shunts)

  println(io, "\n", totalLosses)

  @printf(io, "==========================================================================================================================================================================\n")
  @printf(io, "| %-5s | %-20s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-5s |\n", "Nr", "Bus", "Vn [kV]", "V [kV]", "V [pu]", "phi [deg]", "Pg [MW]", "Qg [MVar]", "Pl [MW]", "Ql [MVar]", "Ps [MW]", "Qs [MVar]", "Type")
  @printf(io, "==========================================================================================================================================================================\n")
  
  pGS = qGS = pLS = qLS = ""
  tpGS = tqGS = tpLS = tqLS = 0.0
  pShunt_str = qShunt_str = ""
  tpShunt = tqShunt = 0.0
  for n in nodes    
    if !isnothing(n._pƩGen) abs(n._pƩGen) > 1e-6
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
    v = n.comp.cVN*n._vm_pu
    nodeName = n.comp.cName
    if !isnothing(n._vmin_pu) && !isnothing(n._vmax_pu)
      if n._vm_pu < n._vmin_pu || n._vm_pu > n._vmax_pu
        nodeName *= " !"
      end
    end

    @printf(io, "| %-5d | %-20s | %-10.1f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-5s |\n", n.busIdx, nodeName, n.comp.cVN, v, n._vm_pu, n._va_deg, pGS, qGS, pLS, qLS, pShunt_str, qShunt_str, typeStr)
  end

  @printf(io, "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  #@show tpGS, tqGS, tpLS, tqLS, tpShunt, tqShunt
  println(io, flowResults)

  if toFile
    close(io)
    #println("Results have been written to $(joinpath(path, filename))")
  end
end


function convertPVtoPQ!(net::ResDataTypes.Net)
  for n in net.nodeVec
    if n._nodeType == ResDataTypes.PV
      busIdx = n.busIdx
      for p in net.prosumpsVec
        if p.comp.cFrom_bus == busIdx          
          setQGenReplacement!(p, n._qƩGen)          
        end
      end
    end
  end
end