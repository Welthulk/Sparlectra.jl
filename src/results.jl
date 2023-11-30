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
  formatted_results  = @sprintf("\n===========================================================================================================================\n")
  formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Branch", "From", "To", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]")
  formatted_results *= @sprintf("===========================================================================================================================\n")

    ∑pv = 0.0
    ∑qv = 0.0
    for br in net.branchVec
        from = br._from
        to = br._to
        bName = br.comp.cName
        pfromVal = (br.fBranchFlow.pFlow === nothing) ? NaN : br.fBranchFlow.pFlow
        qfromVal = (br.fBranchFlow.qFlow === nothing) ? NaN : br.fBranchFlow.qFlow

        ptoVal = (br.tBranchFlow.pFlow === nothing) ? NaN : br.tBranchFlow.pFlow
        qtoVal = (br.tBranchFlow.qFlow === nothing) ? NaN : br.tBranchFlow.qFlow

        VmFrom = (br.fBranchFlow.vm_pu === nothing) ? NaN : br.fBranchFlow.vm_pu
        VaFrom = (br.fBranchFlow.va_deg === nothing) ? NaN : br.fBranchFlow.va_deg

        VmTo = (br.tBranchFlow.vm_pu === nothing) ? NaN : br.tBranchFlow.vm_pu
        VaTo = (br.tBranchFlow.va_deg === nothing) ? NaN : br.tBranchFlow.va_deg

        v1 = VmFrom * exp(im * deg2rad(VaFrom))
        v2 = VmTo * exp(im * deg2rad(VaTo))

        tap = (br.ratio != 0.0) ? br.ratio : 1.0

        Vdiff = v1 / tap - v2

        y = inv((br.r_pu + br.x_pu * im)) #+ 0.5*(br.g_pu + br.b_pu*im)

        Sdiff = Vdiff * conj(Vdiff * y) * net.baseMVA
        pv = real(Sdiff)
        qv = imag(Sdiff)
        ∑pv += abs(pv)
        ∑qv += abs(qv)

        formatted_results *= @sprintf("| %-25s | %-5s | %-5s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f |  %-10.3f|\n", bName, from, to, pfromVal, qfromVal, ptoVal, qtoVal, abs(pv), abs(qv))
    end
    formatted_results *= @sprintf("---------------------------------------------------------------------------------------------------------------------------\n")
    total_losses = @sprintf("total losses (I^2*Z): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

    return formatted_results, total_losses
end


function printACPFlowResults(net::ResDataTypes.Net, ct::Float64, ite::Int, toFile::Bool = false, path::String = "")
  if toFile
    filename = "result_$(net.name).txt"
    file = joinpath(path, filename)
    io = open(file, "w")
    redirect_stdout(io)
  end

  vers = ResDataTypes.SparlectraVersion
  current_date = Dates.format(Dates.now(), "d-u-yy H:M:S")

  formatted_version = format_version(vers)
  flowResults, totalLosses = formatBranchResults(net)

  println("================================================================================")
  println("| SPARLECTRA Version $formatted_version - AC Power Flow Results                          |")
  println("================================================================================")

  vers = ResDataTypes.SparlectraVersion
  current_date = Dates.format(Dates.now(), "d-u-yy H:M:S")

  busses = length(net.nodeVec)
  branches = length(net.branchVec)
  lines = length(net.linesAC)
  trafos = length(net.trafos)
  gens = 0
  loads = 0
  shunts = 0
  auxb = 0

  nodes = sort(net.nodeVec, by = x -> x._kidx)
  npv = 0
  npq = 0
  for n in nodes
    npv += n._nodeType == ResDataTypes.PV ? 1 : 0
    npq += n._nodeType == ResDataTypes.PQ ? 1 : 0
    if occursin("aux_#", n.comp.cName)
      auxb += 1
    end
  end
  for ps in net.prosumpsVec
    loads += ps.proSumptionType == ResDataTypes.Consumption ? 1 : 0
    gens += ps.proSumptionType == ResDataTypes.Injection ? 1 : 0
  end
  shunts = length(net.shuntVec)
  @printf("Date         :%20s\n", current_date)
  @printf("Iterations   :%10d\n", ite)
  @printf("Converged in :%10f seconds\n", ct)
  @printf("Case         :%15s\n", net.name)
  @printf("BaseMVA      :%10d\n", net.baseMVA)
  if auxb > 0
    @printf("Nodes        :%10d (PV: %d PQ: %d (Aux: %d) Slack: %d\n", busses, npv, npq, auxb, 1)
  else
    @printf("Nodes        :%10d (PV: %d PQ: %d Slack: %d)\n", busses, npv, npq, 1)
  end
  @printf("Branches     :%10d\n", branches)
  @printf("Lines        :%10d\n", lines)
  @printf("Trafos       :%10d\n", trafos)
  @printf("Generators   :%10d\n", gens)
  @printf("Loads        :%10d\n", loads)
  @printf("Shunts       :%10d\n", shunts)
  
  println("\n",totalLosses)
  
  @printf("======================================================================================================\n")
  @printf("| %-20s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Bus", "Vn [kV]", "V [pu]", "phi [deg]", "P [MW]", "Q [MVar]", "Type")
  @printf("======================================================================================================\n")
  ∑pLoad = 0.0
  ∑pGen = 0.0
  ∑qLoad = 0.0
  ∑qGen = 0.0
  pVal = 0.0
  qVal = 0.0

  for n in nodes
    if !isnothing(n._pƩGen) && n._pƩGen > 0.0
      pVal = (n._pƩGen === nothing) ? 0.0 : n._pƩGen
      qVal = (n._qƩGen === nothing) ? 0.0 : n._qƩGen
      ∑pGen += pVal
      ∑qGen += qVal
    elseif !isnothing(n._pƩLoad) && n._pƩLoad > 0.0
      pVal = (-n._pƩLoad === nothing) ? 0.0 : -n._pƩLoad
      qVal = (-n._qƩLoad === nothing) ? 0.0 : -n._qƩLoad
      ∑pLoad += pVal
      ∑qLoad += qVal
    end
    typeStr = ResDataTypes.toString(n._nodeType)
    @printf("| %-20s | %-10d | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10s |\n", n.comp.cName, n.comp.cVN, n._vm_pu, n._va_deg, pVal, qVal, typeStr)
  end
  @printf("------------------------------------------------------------------------------------------------------\n")
  println(flowResults)
  if toFile
    # Close the file after redirecting the output
    redirect_stdout(Base.stdout)
    close(io)
    println("Results have been written to $filename")
  end
end
