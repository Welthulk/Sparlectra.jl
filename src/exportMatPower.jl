# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.6.2023
# include-file SparqlQueryCGMES.jl
# further details see: https://github.com/Welthulk/CIMDraw/blob/master/js/matpower.js

const comment = "% This file was generated by Sparlectra exportMatPower.jl\n"


function writeMatpowerCasefile(net::ResDataTypes.Net, filename::String, testcase::String)
  @info "convertion to Matpower CASE-Files, Testcase: ($testcase), Filename: ($filename)"
  base, ext = splitext(filename)
  base = base * ".m"
  NodeDict = Dict{String,ResDataTypes.Node}()
  for n in net.nodeVec
    NodeDict[n.comp.cID] = n
  end

  file = open(filename, "w")
  writeHeader(net.baseMVA, file, testcase)
  writeBusData(net.nodeVec, file)
  WriteGeneratorData(net.baseMVA, NodeDict, net.prosumpsVec, file)

  writeBranchData(net.baseMVA, net.branchVec, file)
  #writeCostData(file)
  close(file)
end # writeMatpowerCasefile

function writeMatpowerCasefile(NodeDict::Dict{String,ResDataTypes.Node}, NodeVec::Vector{ResDataTypes.Node}, ProSumVec::Vector{ResDataTypes.ProSumer}, branchVec::Vector{ResDataTypes.Branch}, testcase::String, sb_mva::Float64, filename::AbstractString, filepath::AbstractString)
  @info "convertion CGMES to Matpower CASE-Files, Testcase: ($testcase)"
  base, ext = splitext(filename)
  base = base * ".m"
  file = open(joinpath(filepath, base), "w")
  writeHeader(sb_mva, file, testcase)
  writeBusData(NodeVec, file)
  WriteGeneratorData(sb_mva, NodeDict, ProSumVec, file)

  writeBranchData(sb_mva, branchVec, file)
  #writeCostData(file)
  close(file)
end # writeMatpowerCasefile

function writeHeader(sb_mva::Float64, file, case::AbstractString)
  write(file, comment)
  write(file, "function mpc = $case\n")
  write(file, "%% MATPOWER Case Format : Version 2\n")
  write(file, "mpc.version = '2';\n")
  write(file, "%%-----  Power Flow Data  -----%%\n")
  write(file, "%% system MVA base\n")
  write(file, "mpc.baseMVA = $sb_mva;\n")
end # writeHeader

function writeBusData(NodeVec::Vector{ResDataTypes.Node}, file)
  write(file, "%% bus data\n")
  write(file, "mpc.bus = [\n")
  write(file, "%bus\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin\n")

  rows = 0
  for node in NodeVec
    rows += 1
    bus_i = node.busIdx
    Name = node.comp.cName
    ID = node.comp.cID
    kIdx = string(node._kidx)

    val = Int(node._nodeType)
    if val <= 0 || val > 4
      @warn "Node $(ID) Node-Type $(val) not supported yet! -> set to 1 (PQ)"
      num = 1
    elseif val == 1
      num = 1
    elseif val == 2
      num = 2
    elseif val == 3
      num = 3
      Name = Name * " (slack)"
    elseif val == 4
      num = 4
      Name = Name * " (isolated)"
    end

    pSum = (node._pƩLoad === nothing) ? 0.0 : node._pƩLoad
    qSum = (node._qƩLoad === nothing) ? 0.0 : node._qƩLoad

    GS = (node._pShunt === nothing) ? 0.0 : node._pShunt
    BS = (node._qShunt === nothing) ? 0.0 : node._qShunt

    area = (node._area === nothing) ? 1 : node._area #https://zepben.github.io/evolve/docs/cim/cim100/TC57CIM/IEC61970/Base/ControlArea/ControlArea/      

    vm_pu = (node._vm_pu === nothing) ? 1.0 : node._vm_pu
    vm_pu = round(vm_pu, digits = 2)
    va_deg = (node._va_deg === nothing) ? 0.0 : node._va_deg
    va_deg = round(va_deg, digits = 2)
    if isnothing(node.comp.cVN)
      @warn "Base Voltage not defined for node $(node.comp.cID) -> set to 1.0"
      baseKV = 1.0
    else
      baseKV = node.comp.cVN
    end

    lZone = (node._lZone === nothing) ? 1 : node._lZone

    Vmax = 1.1
    Vmin = 0.9

    line = string(bus_i, "\t", num, "\t", pSum, "\t", qSum, "\t", GS, "\t", BS, "\t", area, "\t", vm_pu, "\t", va_deg, "\t", baseKV, "\t", lZone, "\t", Vmax, "\t", Vmin, "; %Bus: " * Name * " (ID: " * ID * " kIdx: " * kIdx * ")\n")
    write(file, line)
  end
  write(file, "];\n")
end # writeBusData

function WriteGeneratorData(sb_mva::Float64, NodeDict::Dict{String,ResDataTypes.Node}, ProSumVec::Vector{ResDataTypes.ProSumer}, file)
  write(file, "%% generator data\n")
  write(file, "mpc.gen = [\n")
  write(file, "%bus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin\tPc1\tPc2\tQc1min\tQc1max\tQc2min\tQc2max\tramp_agc\tramp_10\tramp_30\tramp_q\tapf\n")

  for prosum in ProSumVec
    # Motor?
    if isGenerator(prosum.comp) || isExternalNetworkInjection(prosum.comp)
      nodeID = prosum.nodeID
      name = prosum.comp.cName
      cTypStr = string(prosum.comp.cTyp)
      try
        node = NodeDict[nodeID]

        if node.busIdx == 0
          @warn "Node $nodeID has no kidx: $(node.busIdx), $(node)"
          continue
        end

        bus_i = node.busIdx
        Vg = (node._vm_pu === nothing) ? 1.0 : node._vm_pu
        Vg = round(Vg, digits = 2)

        p = (prosum.pVal === nothing) ? 0 : prosum.pVal
        q = (prosum.qVal === nothing) ? 0 : prosum.qVal

        qmax = (prosum.maxQ === nothing) ? 0 : prosum.maxQ
        qmin = (prosum.minQ === nothing) ? 0 : prosum.minQ
        mbase = (prosum.ratedS === nothing) ? sb_mva : prosum.ratedS
        status = 1
        pmax = (prosum.maxP === nothing) ? 0 : prosum.maxP
        pmin = (prosum.minP === nothing) ? 0 : prosum.minP
        pc1 = 0.0
        pc2 = 0.0
        qc1min = 0.0
        qc1max = 0.0
        qc2min = 0.0
        qc2max = 0.0
        ramp_agc = 0.0
        ramp_10 = 0.0
        ramp_30 = 0.0
        ramp_q = 0.0
        apf = 0.0

        line = string(
          bus_i,
          "\t",
          p,
          "\t",
          q,
          "\t",
          qmax,
          "\t",
          qmin,
          "\t",
          Vg,
          "\t",
          mbase,
          "\t",
          status,
          "\t",
          pmax,
          "\t",
          pmin,
          "\t",
          pc1,
          "\t",
          pc2,
          "\t",
          qc1min,
          "\t",
          qc1max,
          "\t",
          qc2min,
          "\t",
          qc2max,
          "\t",
          ramp_agc,
          "\t",
          ramp_10,
          "\t",
          ramp_30,
          "\t",
          ramp_q,
          "\t",
          apf,
          "; %Gen: " * cTypStr * "_" * name * " (ID: " * nodeID * ")\n",
        )
        write(file, line)
      catch err
        @warn "could not handle $nodeID", err
        continue
      end # try
    end # if   
  end # for

  write(file, "];\n")
end # WriteGeneratorData


function writeBranchData(sb_mva::Float64, branchVec::Vector{ResDataTypes.Branch}, file)
  write(file, "%% branch data\n")
  write(file, "mpc.branch = [\n")
  write(file, "%fbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax\n")

  for br in branchVec
    fbus = br.fromBus
    tbus = br.toBus
    r = round(br.r_pu, digits = 8)
    x = round(br.x_pu, digits = 8)
    b = round(br.b_pu, digits = 8)

    rateA = 0 # 0 for unlimited
    rateB = 0 # 0 for unlimited
    rateC = 0 # 0 for unlimited

    ratio = br.ratio
    angle = br.angle
    status = br.status
    angmin = -360 #br.angmin
    angmax = 360 #br.angmax

    type = ResDataTypes.toString(br.comp)
    if br.status == 0 && br.isParallel
      continue
      type = "(parallel & no service) " * type
    elseif br.status == 0
      type = "(no service) " * type
    elseif br.isParallel
      continue
      type = "(parallel) " * type
    end

    cID = br.comp.cID
    line = string(fbus, "\t", tbus, "\t", r, "\t", x, "\t", b, "\t", rateA, "\t", rateB, "\t", rateC, "\t", ratio, "\t", angle, "\t", status, "\t", angmin, "\t", angmax, "; %Branch: " * type * " (ID: " * cID * ")\n")
    write(file, line)
  end # for

  write(file, "];\n")

end # writeBranchData