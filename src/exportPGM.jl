# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 8.2.24
# include-file exportPGM.jl

function exportPGM(net::ResDataTypes.Net, filename::String)
    full_path = filename * ".json"
    @info "export to PGM-Files, Filename: ($full_path)"

    id_counter = 0
    function get_next_id()
        id_counter += 1
        return id_counter
    end

    function get_pgm_nodes(o::Node)
        parameters = OrderedDict{String, Any}()
        @assert isa(o.comp, ImpPGMComp)
        # busId should be origID
        parameters["id"] = get_next_id()
        parameters["u_rated"] = o.comp.cVN * 1e3
        return parameters
    end

    function get_lines(o::ACLineSegment)
        parameters = OrderedDict{String, Any}()
        @assert isa(o.comp, ImpPGMComp)
        parameters["id"] = get_next_id()
        parameters["from_node"] = o.comp.cFrom_bus
        parameters["to_node"] = o.comp.cTo_bus
        parameters["from_status"] = 1
        parameters["to_status"] = 1
        parameters["r1"] = o.r
        parameters["x1"] = o.x
        parameters["c1"] = isnothing(o.c_nf_per_km) ? 0.0 : o.c_nf_per_km * 1e-9
        parameters["tan1"] = isnothing(o.tanδ) ? 0.0 : o.tanδ
        return parameters
    end

    function get_trafos(o::PowerTransformer)
        parameters = OrderedDict{String, Any}()
        @assert isa(o.comp, ImpPGMComp)
        @assert o.isBiWinder == true "Only BiWinder is supported"
        use_default_params = false
        if isnothing(o.exParms)  
          use_default_params = true
          @warn "No additional parameters found"
        end
        
        tap_side = 0
        tap_pos = 0
        tap_min = 0
        tap_max = 0
        tap_nom = 0
        tap_size = 0
        u1 = 0.0
        u2 = 0.0
        winding = nothing
        if o.HVSideNumber == 1
            u1 = o.side1.Vn*1e3
            u2 = o.side2.Vn*1e3
            winding = o.side1
        else
            u1 = o.side2.Vn*1e3
            u2 = o.side1.Vn*1e3
            winding = o.side2
        end

        if o.isControlled
         if o.HVSideNumber == 1
            taps = o.side1.taps
            tap_side = 0
            vn = o.side1.Vn
         else
            taps = o.side2.taps
            tap_side = 1
            vn = o.side2.Vn
         end
         tap_pos = taps.step         
         tap_min = taps.lowStep
         tap_max = taps.highStep
         tap_nom = taps.neutralStep
         #(vs/vn)*100.0

         tap_size = 1e3*taps.voltageIncrement*vn/100.0
        end

        parameters["id"] = get_next_id()
        parameters["from_node"] = o.comp.cFrom_bus
        parameters["to_node"] = o.comp.cTo_bus
        parameters["from_status"] = 1
        parameters["to_status"] = 1
        parameters["u1"] = u1
        parameters["u2"] = u2
        parameters["sn"], parameters["uk"], parameters["pk"], parameters["i0"], parameters["p0"] = use_default_params ? (0.0, 0.0, 0.0, 0.0, 0.0) : (o.exParms.sn, o.exParms.uk, o.exParms.pk, o.exParms.i0, o.exParms.p0)
        parameters["winding_from"] = 1
        parameters["winding_to"] = 1
        parameters["clock"] = 0  # always zero
        parameters["tap_side"] = tap_side
        parameters["tap_pos"] = tap_pos
        parameters["tap_min"] = tap_min
        parameters["tap_max"] = tap_max
        parameters["tap_nom"] = tap_nom
        parameters["tap_size"] = tap_size

        return parameters
    end

    function get_sym_gens(o::ProSumer)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      if o.isAPUNode
        @warn "PV-Node not supported for PGM, set to PQ-Node"
      end  
            
      parameters["id"] = get_next_id()
      parameters["node"] = o.comp.cFrom_bus
      parameters["status"] = 1
      parameters["type"] = 0
      if isnothing(o.pVal)
        @warn "pVal is not set, set to 0.0"
        parameters["p_specified"] = 0.0
      else
        parameters["p_specified"] = o.pVal*1e6
      end

      if isnothing(o.qVal)
        @warn "qVal is not set, set to 0.0"
        parameters["q_specified"] = 0.0
      else
        parameters["q_specified"] = o.qVal*1e6
      end
      
      return parameters
    end

    function get_sym_loads(o::ProSumer)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = get_next_id()
      parameters["node"] = o.comp.cFrom_bus
      parameters["status"] = 1
      parameters["type"] = 0
      parameters["p_specified"] = o.pVal*1e6
      parameters["q_specified"] = o.qVal*1e6

      return parameters
    end

    function get_source(o::Node)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = get_next_id()
      parameters["node"] = o.busIdx
      parameters["status"] = 1 # default
      parameters["u_ref"] = o._vm_pu
      parameters["sk"] = 10000000000.0 # default
      parameters["rx_ratio"] = 0.1 # default
      parameters["z01_ratio"] = 1.0 # default

      return parameters
    end


    function get_shunts(o::Shunt)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = get_next_id()
      parameters["node"] = o.comp.cFrom_bus
      parameters["status"] = 1
      vn = Float64(o.comp.cVN)
      g1, b1 = calcGB_Shunt(o.p_shunt, o.q_shunt, vn)
      parameters["g1"] = g1
      parameters["b1"] = b1

      return parameters
      
    end


    json_data = JSON.json(
        OrderedDict(
            "version" => "1.0",
            "type" => "input",
            "is_batch" => false,
            "attributes" => OrderedDict(),
            "data" => OrderedDict(
                "node" => [get_pgm_nodes(node) for node in net.nodeVec],
                "line" => [get_lines(ac) for ac in net.linesAC],
                "sym_gen" => [get_sym_gens(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Injection && !isSlack(o)],
                "shunt" => [get_shunts(o) for o in net.shuntVec],
                "transformer" => [get_trafos(t) for t in net.trafos],
                "sym_load" => [get_sym_loads(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Consumption],
                "source" => [get_source(o) for o in net.nodeVec if isSlack(o)]
            )
        ),
        2  # specify the indentation level for pretty printing
    )

    open(full_path, "w") do file
        write(file, json_data)
    end
end
