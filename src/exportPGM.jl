# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 8.2.24
# include-file exportPGM.jl

function exportPGM(;net::ResDataTypes.Net, filename::String, useMVATrafoModell::Bool=true)
  full_path = filename * ".json"
  @info "export to PGM-Files, Filename: ($full_path)"

  id_counter = 0
  function get_next_id()
    id_counter += 1
    return id_counter
  end

  function get_pgm_nodes(o::Node)
    parameters = OrderedDict{String,Any}()
    @assert isa(o.comp, ImpPGMComp)
    # busId should be origID
    if o.comp.cTyp == ResDataTypes.AuxBus && !useMVATrafoModell 
      @info "AuxBus: ", o.comp.cName, " ignored..."
      return parameters
    else
      parameters["id"] = get_next_id()
      parameters["u_rated"] = o.comp.cVN * 1e3
      parameters["_cname"] = o.comp.cName
      parameters["_cID"] = o.comp.cID
    end

    return parameters
  end

  function get_lines(o::ACLineSegment)
    parameters = OrderedDict{String,Any}()
    if isa(o.comp, ImpPGMComp)
      parameters["id"] = get_next_id()
      parameters["from_node"] = o.comp.cFrom_bus
      parameters["to_node"] = o.comp.cTo_bus
      parameters["from_status"] = 1
      parameters["to_status"] = 1
      parameters["r1"] = o.r
      parameters["x1"] = o.x
      parameters["c1"] = isnothing(o.c_nf_per_km) ? 0.0 : o.c_nf_per_km * 1e-9
      parameters["tan1"] = isnothing(o.tanδ) ? 0.0 : o.tanδ
      parameters["_cname"] = o.comp.cName
      parameters["_cID"] = o.comp.cID
    else
      @warn "Line $(o.comp.cName) is not exported to PGM (not in service?)"
    end
    return parameters
  end

  function get_trafos(o::PowerTransformer, baseMVA::Float64)
    parameters = OrderedDict{String,Any}()
    #
    if o.isBiWinder == false 
        @warn "Only BiWinder is supported", o.comp.cName
        @show "trafo: ", o
    end  

    if isa(o.comp, ImpPGMComp)
      use_default_params = false

      if isnothing(o.exParms)
        use_default_params = true
        @info "no transforme model data found, recalculation is to be performed"
      end

      tap_side = 0
      tap_pos = 0
      tap_min = 0
      tap_max = 0
      tap_nom = 0
      tap_size = 0
      u1 = 0.0
      u2 = 0.0
      vn = 0.0
      rk = 0.0
      xk = 0.0
      bm = 0.0
      gm = 0.0
      sn = 0.0
      isPerUnitRXGB = false
      ratedU = nothing
      winding = (o._equiParms == 1) ? o.side1 :
           ((o._equiParms == 2) ? o.side2 :
            nothing)
      #@show "winding: ", winding, o.tapSideNumber, o._equiParms
      if !isnothing(winding)
        vn = winding.Vn
        rk = winding.r
        xk = winding.x
        bm = !isnothing(winding.b) ? winding.b : 0.0
        gm = !isnothing(winding.g) ? winding.g : 0.0
        
        ratedU = !isnothing(winding.ratedU) ? winding.ratedU : vn
        sn = !isnothing(winding.ratedS) ? winding.ratedS : 0.0
        isPerUnitRXGB = isPerUnit_RXGB(winding)
      end
      
      
      if o.HVSideNumber == 1
        u1 = o.side1.Vn * 1e3
        u2 = o.side2.Vn * 1e3
      else
        u1 = o.side2.Vn * 1e3
        u2 = o.side1.Vn * 1e3        
      end
      
      if !isnothing(winding) && !isnothing(winding.taps)
        taps = winding.taps
        tap_pos = taps.step
        tap_min = taps.lowStep
        tap_max = taps.highStep
        tap_nom = taps.neutralStep
        #(vs/vn)*100.0

        tap_size = 1e3 * taps.voltageIncrement * vn / 100.0
      end

      parameters["id"] = get_next_id()
      parameters["from_node"] = o.comp.cFrom_bus
      parameters["to_node"] = o.comp.cTo_bus
      parameters["from_status"] = 1
      parameters["to_status"] = 1
      parameters["u1"] = u1
      parameters["u2"] = u2
      if use_default_params
        if sn == 0.0
          @warn "ratedS of transformer $(o.comp.cName) is set to zero, could not perform transformer model data"
          parameters["sn"] = 0.0
          parameters["uk"] = 0.0
          parameters["pk"] = 0.0
          parameters["i0"] = 0.0
          parameters["p0"] = 0.0
        else
          rated_U = !isnothing(ratedU) ? ratedU : vn

          uk, P_kW, i0, Pfe_kW = recalc_trafo_model_data(baseMVA = baseMVA, Sn_MVA = sn, ratedU_kV = rated_U, r_pu = rk, x_pu = xk, b_pu = abs(bm), isPerUnit=isPerUnitRXGB)
          parameters["sn"] = sn * 1e6
          parameters["uk"] = uk
          parameters["pk"] = P_kW * 1e3
          parameters["i0"] = i0
          parameters["p0"] = Pfe_kW * 1e3
        end
      else
        parameters["sn"] = !isnothing(o.exParms.sn) ? o.exParms.sn : 0.0
        parameters["uk"] = !isnothing(o.exParms.uk) ? o.exParms.uk : 0.0
        parameters["pk"] = !isnothing(o.exParms.pk) ? o.exParms.pk : 0.0
        parameters["i0"] = !isnothing(o.exParms.i0) ? o.exParms.i0 : 0.0
        parameters["p0"] = !isnothing(o.exParms.p0) ? o.exParms.p0 : 0.0
      end

      parameters["winding_from"] = 1
      parameters["winding_to"] = 1
      parameters["clock"] = 0  # always zero
      parameters["tap_side"] = tap_side
      parameters["tap_pos"] = tap_pos
      parameters["tap_min"] = tap_min
      parameters["tap_max"] = tap_max
      parameters["tap_nom"] = tap_nom
      parameters["tap_size"] = tap_size
      parameters["_cname"] = o.comp.cName
      parameters["_cID"] = o.comp.cID
    else
      @warn "Transformer $(o.comp.cName) is not exported to PGM (not in service?)"
    end
    return parameters
  end

  function get_sym_gens(o::ProSumer)
    parameters = OrderedDict{String,Any}()
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
      parameters["p_specified"] = abs(o.pVal) * 1e6
    end

    if isnothing(o.qVal)
      @warn "qVal is not set, set to 0.0"
      parameters["q_specified"] = 0.0
    else
      parameters["q_specified"] = abs(o.qVal) * 1e6
    end
    parameters["_cname"] = o.comp.cName
    parameters["_cID"] = o.comp.cID
 
    return parameters
  end

  function get_sym_loads(o::ProSumer)
    parameters = OrderedDict{String,Any}()
    @assert isa(o.comp, ImpPGMComp)

    parameters["id"] = get_next_id()
    parameters["node"] = o.comp.cFrom_bus
    parameters["status"] = 1
    parameters["type"] = 0
    parameters["p_specified"] = abs(o.pVal) * 1e6
    parameters["q_specified"] = abs(o.qVal) * 1e6

    parameters["_cname"] = o.comp.cName
    parameters["_cID"] = o.comp.cID

    return parameters
  end

  function get_source(o::Node)
    parameters = OrderedDict{String,Any}()
    @assert isa(o.comp, ImpPGMComp)

    parameters["id"] = get_next_id()
    parameters["node"] = o.busIdx
    parameters["status"] = 1 # default
    if isnothing(o._vm_pu)
      @warn "u_ref is not set, set to 1.0"
      parameters["u_ref"] = 1.0
    else
      parameters["u_ref"] = o._vm_pu
    end
    parameters["sk"] = 10000000000.0 # default
    parameters["rx_ratio"] = 0.1 # default
    parameters["z01_ratio"] = 1.0 # default
    parameters["_cname"] = o.comp.cName
    parameters["_cID"] = o.comp.cID

    return parameters
  end

  function get_shunts(o::Shunt)
    parameters = OrderedDict{String,Any}()
    @assert isa(o.comp, ImpPGMComp)

    parameters["id"] = get_next_id()
    parameters["node"] = o.comp.cFrom_bus
    parameters["status"] = 1
    vn = Float64(o.comp.cVN)
    g1, b1 = calcGB_Shunt(o.p_shunt, o.q_shunt, vn)
    parameters["g1"] = g1
    parameters["b1"] = b1
    parameters["_cname"] = o.comp.cName
    parameters["_cID"] = o.comp.cID


    return parameters
  end

  json_data = JSON.json(
    OrderedDict(
      "version" => "1.0",
      "type" => "input",
      "is_batch" => false,
      "attributes" => OrderedDict(),
      "data" => OrderedDict(
        "node" => filter(x -> !isempty(x), [get_pgm_nodes(node) for node in net.nodeVec]),
        "line" => filter(x -> !isempty(x), [get_lines(ac) for ac in net.linesAC]),
        "sym_gen" => filter(x -> !isempty(x), [get_sym_gens(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Injection && !isSlack(o)]),
        "shunt" => filter(x -> !isempty(x), [get_shunts(o) for o in net.shuntVec]),
        "transformer" => filter(x -> !isempty(x), [get_trafos(t, net.baseMVA) for t in net.trafos]),
        "sym_load" => filter(x -> !isempty(x), [get_sym_loads(o) for o in net.prosumpsVec if o.proSumptionType == ResDataTypes.Consumption]),
        "source" => filter(x -> !isempty(x), [get_source(o) for o in net.nodeVec if isSlack(o)]),
      ),
    ),
    2,  # specify the indentation level for pretty printing
  )

  open(full_path, "w") do file
    write(file, json_data)
  end
end
