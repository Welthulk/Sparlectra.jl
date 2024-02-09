# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 8.2.24
# include-file exportPGM.jl

function exportPGM(net::ResDataTypes.Net, filename::String)
    full_path = filename * ".json"
    @info "export to PGM-Files, Filename: ($full_path)"

    function get_pgm_node_dict(comp::ImpPGMComp)
        parameters = OrderedDict{String, Any}()
        @assert isa(comp, ImpPGMComp)
        parameters["id"] = comp.cOrigId
        parameters["u_rated"] = comp.cVN * 1e3
        return parameters
    end

    function get_line_parameters(ac::ACLineSegment)
        parameters = OrderedDict{String, Any}()
        @assert isa(ac.comp, ImpPGMComp)
        parameters["id"] = ac.comp.cOrigId
        parameters["from_node"] = ac.comp.cFrom_bus
        parameters["to_node"] = ac.comp.cTo_bus
        parameters["from_status"] = 1
        parameters["to_status"] = 1
        parameters["r1"] = ac.r
        parameters["x1"] = ac.x
        parameters["c1"] = isnothing(ac.c_nf_per_km) ? 0.0 : ac.c_nf_per_km * 1e-9
        parameters["tan1"] = isnothing(ac.tanδ) ? 0.0 : ac.tanδ
        return parameters
    end

    function get_trafo_parameters(t::PowerTransformer, baseMVA::Float64)
        parameters = OrderedDict{String, Any}()
        @assert isa(t.comp, ImpPGMComp)
        @assert t.isBiWinder == true "Only BiWinder is supported"
        @assert !isnothing(t.exParms)  "No additional parameters found"
        
        tap_side = 0
        tap_pos = 0
        tap_min = 0
        tap_max = 0
        tap_nom = 0
        tap_size = 0
        u1 = 0.0
        u2 = 0.0
        winding = nothing
        if t.HVSideNumber == 1
            u1 = t.side1.Vn*1e3
            u2 = t.side2.Vn*1e3
            winding = t.side1
        else
            u1 = t.side2.Vn*1e3
            u2 = t.side1.Vn*1e3
            winding = t.side2
        end

        if t.isControlled
         if t.HVSideNumber == 1
            taps = t.side1.taps
            tap_side = 0
            vn = t.side1.Vn
         else
            taps = t.side2.taps
            tap_side = 1
            vn = t.side2.Vn
         end
         tap_pos = taps.step         
         tap_min = taps.lowStep
         tap_max = taps.highStep
         tap_nom = taps.neutralStep
         #(vs/vn)*100.0

         tap_size = 1e3*taps.voltageIncrement*vn/100.0
        end

        parameters["id"] = t.comp.cOrigId
        parameters["from_node"] = t.comp.cFrom_bus
        parameters["to_node"] = t.comp.cTo_bus
        parameters["from_status"] = 1
        parameters["to_status"] = 1
        parameters["u1"] = u1
        parameters["u2"] = u2
        parameters["sn"] = t.exParms.sn
        parameters["uk"] = t.exParms.uk
        parameters["pk"] = t.exParms.pk
        parameters["i0"] = t.exParms.i0
        parameters["p0"] = t.exParms.p0
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

    function get_sym_gen(o::ProSumer)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = o.comp.cOrigId
      parameters["node"] = o.comp.cFrom_bus
      parameters["status"] = 1
      parameters["type"] = 0
      parameters["p_specified"] = o.pVal*1e6
      parameters["q_specified"] = o.qVal*1e6
      return parameters
    end

    function get_sym_load(o::ProSumer)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = o.comp.cOrigId
      parameters["node"] = o.comp.cFrom_bus
      parameters["status"] = 1
      parameters["type"] = 0
      parameters["p_specified"] = o.pVal*1e6
      parameters["q_specified"] = o.qVal*1e6

      return parameters
    end
    function get_source(o::ProSumer)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = o.comp.cOrigId
      parameters["node"] = o.referencePri
      parameters["status"] = 1 # default
      parameters["u_ref"] = o.vm_pu
      parameters["sk"] = 10000000000.0 # default
      parameters["rx_ratio"] = 0.1 # default
      parameters["z01_ratio"] = 1.0 # default

      return parameters
    end


    function get_shunt(o::Shunt)
      parameters = OrderedDict{String, Any}()
      @assert isa(o.comp, ImpPGMComp)

      parameters["id"] = o.comp.cOrigId
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
                "node" => [get_pgm_node_dict(node.comp) for node in net.nodeVec],
                "line" => [get_line_parameters(ac) for ac in net.linesAC],
                "sym_gen" => [get_sym_gen(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Injection && !isSlack(o)],
                "shunt" => [get_shunt(o) for o in net.shuntVec],
                "transformer" => [get_trafo_parameters(t, net.baseMVA) for t in net.trafos],
                "sym_load" => [get_sym_load(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Consumption],
                "source" => [get_source(o) for o in net.prosumpsVec if isSlack(o)]
            )
        ),
        2  # specify the indentation level for pretty printing
    )

    open(full_path, "w") do file
        write(file, json_data)
    end
end
