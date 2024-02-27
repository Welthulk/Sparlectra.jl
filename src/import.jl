# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 28.6.2023
# include-file import.jl

function jsonparser(filename, debug::Bool = false)
  json_data = read(filename, String)
  data_dict = JSON.parse(json_data)
  netName = data_dict["Net"]
  baseMVA = data_dict["BaseMVA"]

  buses = data_dict["Buses"]
  idx = 0
  for bus in buses
    idx += 1
    if haskey(bus, "bus")
      if bus["bus"] == ""
        msg = "bus number is empty!"
        throw(msg)
      else
        idx = bus["bus"]
      end
    else
      msg = "No bus number in this net!"
      throw(msg)
    end

    if bus["name"] == ""
      bus["name"] = "Bus_#" * string(idx)
    end

    if bus["id"] == ""
      bus["id"] = UUIDs.uuid4()
    end
  end

  linemod = nothing
  try
    linemod = data_dict["LineModelling"]
    #idx=0
    for lm in linemod
      idx += 1
      #if line["model"]== ""            
      #    line["model"] = "LineMod_#" * string(idx)
      #end
      if lm["id"] == ""
        lm["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No LineModelling in this net")
    end
  end

  lines = data_dict["ACLines"]
  #idx=0
  for line in lines
    idx += 1
    #if line["name"]== ""            
    #    line["name"] = "ACL_#" * string(idx)
    #end
    if line["id"] == ""
      line["id"] = UUIDs.uuid4()
    end
  end

  vkDepChr = nothing
  try
    vkDepChr = data_dict["VK-Characteristics"]
  catch
    if debug
      println("No VK-Characteristics in this net")
    end
  end

  tapMod = nothing
  try
    tapMod = data_dict["TapChangerModelling"]
    idx = 0
    for tap in tapMod
      idx += 1
      if tap["model"] == ""
        tap["model"] = "TapMod_#" * string(idx)
      end
      if tap["id"] == ""
        tap["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No TapChangerModelling in this net")
    end
  end

  trafoWT2 = nothing
  try
    trafoWT2 = data_dict["TwoWindingTransformers"]
    idx = 0
    for trafo in trafoWT2
      idx += 1
      if trafo["name"] == ""
        trafo["name"] = "2WT_#" * string(idx)
      end
      if trafo["id"] == ""
        trafo["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No TwoWindingTransformers in this net")
    end
  end

  trafoWT3 = nothing
  try
    trafoWT3 = data_dict["ThreeWindingTransformers"]
    idx = 0
    for trafo in trafoWT3
      idx += 1
      if trafo["name"] == ""
        trafo["name"] = "3WT_#" * string(idx)
      end
      if trafo["id"] == ""
        trafo["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No ThreeWindingTransformers in this net")
    end
  end

  shunts = nothing
  try
    shunts = data_dict["Shunts"]
    idx = 0
    for shunt in shunts
      idx += 1
      if shunt["name"] == ""
        shunt["name"] = "SH_#" * string(idx)
      end
      if shunt["id"] == ""
        shunt["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No Shunts in this net")
    end
  end

  loads = nothing
  try
    loads = data_dict["Loads"]
    idx = 0
    for load in loads
      idx += 1
      if load["name"] == ""
        load["name"] = "LD_#" * string(idx)
      end
      if load["id"] == ""
        load["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No Loads in this net")
    end
  end
  sgens = nothing
  try
    sgens = data_dict["StaticGenerators"]
    idx = 0
    for sgen in sgens
      idx += 1
      if sgen["name"] == ""
        sgen["name"] = "SGEN_#" * string(idx)
      end
      if sgen["id"] == ""
        sgen["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No StaticGenerators in this net")
    end
  end
  vgens = nothing
  try
    vgens = data_dict["VoltageControlledGenerators"]
    idx = 0
    for vgen in vgens
      idx += 1
      if vgen["name"] == ""
        vgen["name"] = "VGEN_#" * string(idx)
      end
      if vgen["id"] == ""
        vgen["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No VoltageControlledGenerators in this net")
    end
  end
  ext_grids = nothing
  try
    ext_grids = data_dict["ExternalGrids"]
    idx = 0
    for ext_grid in ext_grids
      idx += 1
      if ext_grid["name"] == ""
        ext_grid["name"] = "ExtGrid_#" * string(idx)
      end
      if ext_grid["id"] == ""
        ext_grid["id"] = UUIDs.uuid4()
      end
    end
  catch
    if debug
      println("No ExternalGrids in this net")
    end
  end

  if debug
    println("Net: ", netName, "\n")
    println("BaseMVA: ", baseMVA, "\n")
    println("Buses: ", buses, "\n")
    println("linemod: ", linemod, "\n")
    println("Lines: ", lines, "\n")
    println("VkDep: ", vkDepChr, "\n")
    println("TapMod: ", tapMod, "\n")
    println("Trafowt2: ", trafoWT2, "\n")
    println("Trafowt3: ", trafoWT3, "\n")
    println("Shunts: ", shunts, "\n")
    println("Loads: ", loads, "\n")
    println("Sgens: ", sgens, "\n")
    println("Vgens: ", vgens, "\n")
    println("Ext Grids: ", ext_grids, "\n")
  end

  return netName, baseMVA, Dict("buses" => buses, "linemod" => linemod, "lines" => lines, "vkDepChr" => vkDepChr, "tapmod" => tapMod, "trafowt2" => trafoWT2, "trafowt3" => trafoWT3, "shunts" => shunts, "loads" => loads, "sgens" => sgens, "vgens" => vgens, "ext_grids" => ext_grids)
end

# brute force reading of matpower casefiles
function casefileparser(filename)
  function processLine(line::AbstractString)::Vector{Float64}
    line = chop(line)  # Run chop here to remove unnecessary characters at the end of the line
    if endswith(line, ';')
      line = chop(line[1:end-1])  # Remove semicolon at end if present
    end
    try

      # Split by tabs and spaces
      values = split(line, ['\t', ' '])

      # Removing empty strings
      values = filter(x -> x ≠ "", values)

      # Parsing the values ​​as Float64
      values = parse.(Float64, values)

      if isempty(values)
        println("Error in line: ", line)
        return []
      end
      return values
    catch err
      println("Error: ", err)
      println("Error in line: ", line)
      return []
    end
  end

  file = open(filename, "r")

  bus_data_block = ""
  gen_data_block = ""
  branch_data_block = ""
  len_bus = 0
  len_gen = 0
  len_branch = 0
  mpc_version = 0
  case_name = ""
  baseMVA = 0.0
  block = 1
  for line in eachline(file)
    if occursin("function mpc = ", line)
      case_name = line
      case_name = split(case_name, "=")
      case_name = case_name[2]
    elseif occursin("mpc.version", line)
      version = chop(line)
      version = split(version, "'")
      mpc_version = parse(Int, version[2])
      mpc_version != 2 && error("MPC version $(mpc_version) is not supported")
    elseif occursin("mpc.baseMVA", line)
      baseMVA = chop(line)
      baseMVA = split(baseMVA, "=")
      baseMVA = parse(Float64, baseMVA[2])
    elseif occursin("mpc.bus = [", line) && block == 1
      bus_data_block = line
      block = 1
    elseif occursin("];", line)
      if block == 1
        bus_data_block *= "\n$line"
        block = 2
      elseif block == 2
        gen_data_block *= "\n$line"
        block = 3
      elseif block == 3
        branch_data_block *= "\n$line"
        break
      end
    elseif occursin("%", line)
      continue
    elseif !isempty(bus_data_block) && block == 1
      # Add the line to the block
      bus_data_block *= "\n$line"
      len_bus += 1
    elseif !isempty(gen_data_block) && block == 2
      # Add the line to the block
      gen_data_block *= "\n$line"
      len_gen += 1
    elseif !isempty(branch_data_block) && block == 3
      # Add the line to the block
      branch_data_block *= "\n$line"
      len_branch += 1
    elseif occursin("mpc.gen = [", line) && block == 2
      # Start collecting the block from this timele            
      gen_data_block = line
    elseif occursin("mpc.branch = [", line) && block == 3
      branch_data_block = line
    end
  end

  close(file)

  isempty(bus_data_block) && error("Pattern not found")

  bus_data_block = replace(bus_data_block, r"#.*\n" => "\n")

  bus_zeilen = split(bus_data_block, '\n')
  gen_zeilen = split(gen_data_block, '\n')
  branch_zeilen = split(branch_data_block, '\n')

  mpc_bus = zeros(Float64, len_bus, 13)
  mpc_gen = zeros(Float64, len_gen, 21)
  mpc_branch = zeros(Float64, len_branch, 13)

  anz = 0
  for (i, zeile) in enumerate(bus_zeilen)
    if occursin("[", zeile)
      continue
    elseif occursin("];", zeile)
      break
    else
      werte = processLine(zeile)
      mpc_bus[i-1, :] = werte[1:13]
      anz += 1
    end#
  end

  for (i, zeile) in enumerate(gen_zeilen)
    if occursin("[", zeile)
      continue
    elseif occursin("];", zeile)
      break
    else
      werte = processLine(zeile)
      mpc_gen[i-1, :] = werte[1:21]
    end#
  end

  for (i, zeile) in enumerate(branch_zeilen)
    if occursin("[", zeile)
      continue
    elseif occursin("];", zeile)
      break
    else
      werte = processLine(zeile)
      mpc_branch[i-1, :] = werte[1:13]
    end
  end

  # sorting the bus matrix
  # extract the bus number columne
  busidx = mpc_bus[:, 1]

  # get the original bus number orders
  busSequence = sortperm(busidx)

  # apply the order to the mpc_bus array
  mpc_bus = mpc_bus[busSequence, :]
  
  return case_name, baseMVA, mpc_bus, mpc_gen, mpc_branch
end

function pgmparser(filename)
  json_data = read(filename, String)
  data_dict = JSON.parse(json_data)
  
  version = data_dict["version"]
  if version != "1.0"
    msg = "Version $(version) of PGM-File not supported!"
    throw(msg)
  end
  
  nodes = data_dict["data"]["node"]
  lines = data_dict["data"]["line"]
  if haskey(data_dict["data"], "transformer")
    wt2 = data_dict["data"]["transformer"]
  else
    wt2 = []
  end
   
  wt3 = nothing
  sym_gens = data_dict["data"]["sym_gen"]
  if haskey(data_dict["data"], "sym_load")
    sym_loads = data_dict["data"]["sym_load"]  
  else
    sym_loads = []
  end    
    
  if haskey(data_dict["data"], "shunt")
    shunts = data_dict["data"]["shunt"]
  else
    shunts = []
  end  
  source = data_dict["data"]["source"]
  return nodes, lines, wt2, wt3, sym_gens, sym_loads, shunts, source
  
end
