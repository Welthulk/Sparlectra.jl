# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 8.2.24
# include-file exportPGM.jl


function exportPGM(net::ResDataTypes.Net, filename::String)
  full_path = filename*".json"
  @info "export to PGM-Files, Filename: ($full_path)"

  
  function get_pgm_node_dict(comp::ImpPGMComp)
    parameters = Dict{String, Any}()
    @assert isa(comp, ImpPGMComp)
    parameters["id"] = comp.cOrigId
    parameters["u_rated"] = comp.cVN*1e3    
    return parameters
  end
  # Funktion, um die Line-Parameter im gewünschten Format zu generieren
  function get_line_parameters(ac::ACLineSegment)
    parameters = Dict{String, Any}()
    @assert isa(ac.comp, ImpPGMComp)
    parameters["id"] = ac.comp.cOrigId
    parameters["from_node"] = ac.comp.cFrom_bus
    parameters["to_node"] = ac.comp.cTo_bus
    parameters["from_status"] = 1
    parameters["to_status"] = 1
    parameters["r1"] = ac.r
    parameters["x1"] = ac.x
    parameters["c1"] = isnothing(ac.c_nf_per_km) ? 0.0 : ac.c_nf_per_km*1e-9
    parameters["tan1"] = isnothing(ac.tanδ) ? 0.0 : ac.tanδ
    #parameters["r0"] = 0.0  # not used
    #parameters["x0"] = 0.0  # not used
    #parameters["c0"] = 0.0  # not used
    #parameters["tan0"] = 0.0  # not used
    #parameters["i_n"] = 10000.0  # not used

    return parameters
  end

  function get_trafo_parameters(t::PowerTransformer, baseMVA::Float64)
    parameters = Dict{String, Any}()
    @assert isa(t.comp, ImpPGMComp)
    @assert t.isBiWinder == true "Only BiWinder is supported"
    tap_side = 0
    tap_pos = 0
    tap_min = 0
    tap_max = 0
    tap_nom = 0
    tap_size = 0
    u1 = 0.0
    u2 = 0.0
    sn = baseMVA*1e9
    uk = 0.0
    pk = 0.0
    i0 = 0.0
    p0 = 0.0
    if t.isControlled
      if t.HVSideNumber == 1
        taps = t.side1.taps

      else
        taps = t.side2.taps
        tap_side = 1
      end
    end
    
    parameters["id"] = t.comp.cOrigId
    parameters["from_node"] = t.comp.cFrom_bus
    parameters["to_node"] = t.comp.cTo_bus
    parameters["from_status"] = 1
    parameters["to_status"] = 1    
    parameters["u1"] = u1
    parameters["u2"] = u2    
    parameters["sn"] = sn    
    parameters["uk"] = uk
    parameters["pk"] = pk
    parameters["i0"] = i0
    parameters["p0"] = p0    
    parameters["winding_from"] = 1
    parameters["winding_to"] = 1
    parameters["clock"] = 0 # always zero    
    parameters["tap_side"] = tap_side
    parameters["tap_pos"] = tap_pos
    parameters["tap_min"] = tap_min
    parameters["tap_max"] = tap_max
    parameters["tap_nom"] = tap_nom
    parameters["tap_size"] = tap_size


    return parameters
  end

  
  json_data = JSON.json(
    Dict(
        "version" => "1.0",
        "type" => "input",
        "is_batch" => false,
        "attributes" => Dict(),
        "data" => Dict(
            "node" => [get_pgm_node_dict(node.comp) for node in net.nodeVec],            
            "line" => [get_line_parameters(ac) for ac in net.linesAC],
            "transformer": [get_trafo_parameters(t, net.baseMVA) for t in net.trafos]                
            #"sym_gen" => sym_gens,
            #"shunt" => shunts,
            #"transformer" => transformers,
            #"sym_load" => sym_loads,
            #"source" => sources
        )
    )
 )
open(full_path, "w") do file
    write(file, json_data)
end



end # exportPGM