# Copyright 2023–2025 Udo Schmitz
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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 28.6.2023
# include-file import.jl

# Parser for MATPOWER case files
"""
    casefileparser(filename)

Parses a MATPOWER case file and returns the case name, base power, and bus, generator, and branch data.

# Arguments
- `filename`: The name of the file to parse.

# Returns
- `case_name`: The name of the case.
- `baseMVA`: The base power in MVA.
- `mpc_bus`: The bus data matrix.
- `mpc_gen`: The generator data matrix.
- `mpc_branch`: The branch data matrix.

# Example
```julia
case_name, baseMVA, mpc_bus, mpc_gen, mpc_branch = casefileparser("case9.m")
```
"""
function casefileparser(filename)
  function processLine(line::AbstractString)::Vector{Float64}
    line = chop(line)  # Run chop here to remove unnecessary characters at the end of the line
    # Ignore comments after ;
    if endswith(line, ';')
      line = chop(line[1:(end-1)])  # Remove semicolon at end if present
    end

    try
      # Split by tabs and spaces
      values = split(line, ['\t', ' '])

      # Removing empty strings
      values = filter(x -> x ≠ "", values)

      # Parsing the values as Float64
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
    elseif isempty(strip(line))
      continue
    elseif startswith(line, "%")
      continue
    elseif !isempty(bus_data_block) && block == 1
      # Add the line to the block
      bus_data_block *= "\n$line"
      len_bus += 1
    elseif !isempty(gen_data_block) && block == 2
      # Add the line to the block
      gen_data_block *= "\n$line"
      len_gen += 1
    elseif !isempty(strip(branch_data_block)) && block == 3
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

  isempty(bus_data_block) && error("Bus Pattern not found")
  isempty(gen_data_block) && error("Gen Pattern not found")
  isempty(branch_data_block) && error("Branch Pattern not found")

  bus_data_block = replace(bus_data_block, r"%.*\n" => "\n")
  branch_data_block = replace(branch_data_block, r"%.*\n" => "\n")
  gen_data_block = replace(gen_data_block, r"%.*\n" => "\n")

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

  return strip(case_name), baseMVA, mpc_bus, mpc_gen, mpc_branch
end

export casefileparser
