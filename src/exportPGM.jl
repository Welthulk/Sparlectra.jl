# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 8.2.24
# include-file exportPGM.jl

function exportPGM(; net::ResDataTypes.Net, filename::String, useMVATrafoModell::Bool, exportSlackGen::Bool = false)
  full_path = strip(filename) * ".json"
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
      parameters["r1"] = o.r * o.length
      parameters["x1"] = o.x * o.length
      parameters["c1"] = isnothing(o.c_nf_per_km) ? 0.0 : o.c_nf_per_km * 1e-9 * o.length
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

    if o.isBiWinder && isa(o.comp, ImpPGMComp)
      has_model_data = true

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
      winding = (o._equiParms == 1) ? o.side1 : ((o._equiParms == 2) ? o.side2 : nothing)
      @assert o._equiParms > 0
      if !isnothing(winding)
        vn = winding.Vn
        rk = winding.r
        xk = winding.x
        bm = !isnothing(winding.b) ? winding.b : 0.0
        gm = !isnothing(winding.g) ? winding.g : 0.0

        ratedU = !isnothing(winding.ratedU) ? winding.ratedU : vn
        sn = !isnothing(winding.ratedS) ? winding.ratedS : 0.0
        isPerUnitRXGB = isPerUnit_RXGB(winding)
        if isnothing(winding.modelData)
          has_model_data = false
          @info "no transforme model data found, recalculation is to be performed"
        end
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

        tap_size = taps.voltageIncrement_kV * 1e3
      end

      parameters["id"] = get_next_id()
      parameters["from_node"] = o.comp.cFrom_bus
      parameters["to_node"] = o.comp.cTo_bus
      parameters["from_status"] = 1
      parameters["to_status"] = 1
      parameters["u1"] = u1
      parameters["u2"] = u2
      if !has_model_data
        if sn == 0.0
          @warn "ratedS of transformer $(o.comp.cName) is set to zero, could not perform transformer model data"
          parameters["sn"] = 0.0
          parameters["uk"] = 0.0
          parameters["pk"] = 0.0
          parameters["i0"] = 0.0
          parameters["p0"] = 0.0
        else
          rated_U = !isnothing(ratedU) ? ratedU : vn
          uk, P_kW, i0, Pfe_kW = recalc_trafo_model_data(baseMVA = baseMVA, Sn_MVA = sn, ratedU_kV = rated_U, r_pu = rk, x_pu = xk, b_pu = abs(bm), isPerUnit = isPerUnitRXGB)
          parameters["sn"] = sn * 1e6
          parameters["uk"] = uk
          parameters["pk"] = P_kW * 1e3
          parameters["i0"] = i0
          parameters["p0"] = Pfe_kW * 1e3
        end
      else
        parameters["sn"] = !isnothing(winding.modelData.sn_MVA) ? winding.modelData.sn_MVA * 1e6 : 0.0
        parameters["uk"] = !isnothing(winding.modelData.vk_percent) ? winding.modelData.vk_percent * 1e-2 : 0.0
        parameters["pk"] = !isnothing(winding.modelData.pk_kW) ? winding.modelData.pk_kW * 1e3 : 0.0
        parameters["i0"] = !isnothing(winding.modelData.i0_percent) ? winding.modelData.i0_percent * 1e-2 : 0.0
        parameters["p0"] = !isnothing(winding.modelData.p0_kW) ? winding.modelData.p0_kW * 1e3 : 0.0
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

  function get_3WTtrafos(o::PowerTransformer, baseMVA::Float64)
    parameters = OrderedDict{String,Any}()
    if isa(o.comp, ImpPGMComp3WT)

      #@show "trafo: ", o
      u1   = !isnothing(o.side1.ratedU) ? o.side1.ratedU * 1e3 : o.side1.Vn * 1e3
      u2   = !isnothing(o.side2.ratedU) ? o.side2.ratedU * 1e3 : o.side2.Vn * 1e3
      u3   = o.side3.ratedU * 1e3
      sn_1 = !isnothing(o.side1.ratedS) ? o.side1.ratedS * 1e6 : 0.0
      sn_2 = !isnothing(o.side2.ratedS) ? o.side2.ratedS * 1e6 : 0.0
      sn_3 = !isnothing(o.side3.ratedS) ? o.side3.ratedS * 1e6 : 0.0

      uk_12 = 0.0
      uk_13 = 0.0
      uk_23 = 0.0
      pk_12 = 0.0
      pk_13 = 0.0
      pk_23 = 0.0
      i0 = 0.0
      p0 = 0.0

      #@show o.side1.modelData
      #@show o.side2.modelData
      #@show o.side3.modelData

      #@assert !(isnothing(o.side1.modelData) || isnothing(o.side2.modelData) || isnothing(o.side3.modelData)) "modelData not supported for 3WT-trafo"

      l = 0
      for side in [o.side1, o.side2, o.side3]
        l += 1
        if isnothing(side.modelData)
          @info "no transforme model data found, recalculation is to be performed"          
          bm = isnothing(side) ? 0.0 : side.b
          uk, P_kW, _i0, _Pfe_kW = recalc_trafo_model_data(baseMVA = baseMVA, Sn_MVA = side.ratedS, ratedU_kV = side.ratedU, r_pu = side.r, x_pu = side.x, b_pu = abs(bm), isPerUnit = side.isPu_RXGB)

          if _Pfe_kW * 1e3 > p0
            p0 = _Pfe_kW * 1e3
          end

          if _i0 > i0
            i0 = _i0
          end

          if l == 1
            uk_12 = uk
            pk_12 = P_kW * 1e3
          elseif l == 2
            uk_13 = uk
            pk_13 = P_kW * 1e3
          else
            uk_23 = uk
            pk_23 = P_kW * 1e3
          end
        else
          if l == 1
            uk_12 = side.modelData.vk_percent * 1e-2
            pk_12 = side.modelData.pk_kW * 1e3
            _i0 = side.modelData.i0_percent * 1e-2
            _p0 = side.modelData.p0_kW * 1e3
            i0 = max(i0, _i0)
            p0 = max(p0, _p0)
          elseif l == 2
            uk_13 = side.modelData.vk_percent * 1e-2
            pk_13 = side.modelData.pk_kW * 1e3
            _i0 = side.modelData.i0_percent * 1e-2
            _p0 = side.modelData.p0_kW * 1e3
            i0 = max(i0, _i0)
            p0 = max(p0, _p0)
          else
            uk_23 = side.modelData.vk_percent * 1e-2
            pk_23 = side.modelData.pk_kW * 1e3
            _i0 = side.modelData.i0_percent * 1e-2
            _p0 = side.modelData.p0_kW * 1e3
            i0 = max(i0, _i0)
            p0 = max(p0, _p0)
          end
        end
      end

      winding_1 = 0
      winding_2 = 0
      winding_3 = 0
      clock_12  = 0
      clock_13  = 0

      tap_side = 0
      tap_pos = 0
      tap_min = 0
      tap_max = 0
      tap_nom = 0
      tap_size = 0

      if o.tapSideNumber == 1
        tap_side = 0
        tap_pos = o.side1.taps.step
        tap_min = o.side1.taps.lowStep
        tap_max = o.side1.taps.highStep
        tap_nom = o.side1.taps.neutralStep
        tap_size = 1e3 * o.side1.taps.voltageIncrement * o.side1.Vn / 100.0
      elseif o.tapSideNumber == 2
        tap_side = 1
        tap_pos = o.side2.taps.step
        tap_min = o.side2.taps.lowStep
        tap_max = o.side2.taps.highStep
        tap_nom = o.side2.taps.neutralStep
        tap_size = 1e3 * o.side2.taps.voltageIncrement * o.side2.Vn / 100.0
      elseif o.tapSideNumber == 3
        tap_side = 2
        tap_pos = o.side3.taps.step
        tap_min = o.side3.taps.lowStep
        tap_max = o.side3.taps.highStep
        tap_nom = o.side3.taps.neutralStep
        tap_size = 1e3 * o.side3.taps.voltageIncrement * o.side3.Vn / 100.0
      else
        tap_side = 0
        tap_pos = 0
        tap_min = 0
        tap_max = 0
        tap_nom = 0
        tap_size = 0
      end

      node_1 = o.comp.cHV_bus
      node_2 = o.comp.cMV_bus
      node_3 = o.comp.cLV_bus

      status_1 = 1
      status_2 = 1
      status_3 = 1
      parameters["id"] = get_next_id()
      parameters["u1"] = u1
      parameters["u2"] = u2
      parameters["u3"] = u3
      parameters["sn_1"] = sn_1
      parameters["sn_2"] = sn_2
      parameters["sn_3"] = sn_3
      parameters["uk_12"] = uk_12
      parameters["uk_13"] = uk_13
      parameters["uk_23"] = uk_23
      parameters["pk_12"] = pk_12
      parameters["pk_13"] = pk_13
      parameters["pk_23"] = pk_23
      parameters["i0"] = i0
      parameters["p0"] = p0
      parameters["winding_1"] = winding_1
      parameters["winding_2"] = winding_2
      parameters["winding_3"] = winding_3
      parameters["clock_12"] = clock_12
      parameters["clock_13"] = clock_13
      parameters["tap_side"] = tap_side
      parameters["tap_pos"] = tap_pos
      parameters["tap_min"] = tap_min
      parameters["tap_max"] = tap_max
      parameters["tap_nom"] = tap_nom
      parameters["tap_size"] = tap_size
      parameters["node_1"] = node_1
      parameters["node_2"] = node_2
      parameters["node_3"] = node_3
      parameters["status_1"] = status_1
      parameters["status_2"] = status_2
      parameters["status_3"] = status_3
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


    parameters["id"] = get_next_id()
    parameters["node"] = o.comp.cFrom_bus
    parameters["status"] = 1
    parameters["type"] = 0

    if isnothing(o.pVal)
      @warn "pVal is not set, set to 0.0"
      parameters["p_specified"] = 0.0
    else
      parameters["p_specified"] = o.pVal * 1e6
    end

    if isnothing(o.qVal)
      @warn "qVal is not set, set to 0.0"
      parameters["q_specified"] = 0.0
    else
      parameters["q_specified"] = o.qVal * 1e6
    end
    if isAPUNode(o)
      parameters["_c"] = "converted to PQ-Node"
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
    parameters["p_specified"] = o.pVal * 1e6
    parameters["q_specified"] = o.qVal * 1e6

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
    parameters["sk"] = 1e40  # minimize Voltage lost
    parameters["rx_ratio"] = 0.0 # should >= 0.0
    parameters["z01_ratio"] = 0.1 # should > 0.0
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
    g1, b1 = getGBShunt(o)    
    parameters["g1"] = g1
    parameters["b1"] = b1
    parameters["_cname"] = o.comp.cName
    parameters["_cID"] = o.comp.cID

    return parameters
  end

  function filter_gens(o::ProSumer)
    if exportSlackGen
      return o.proSumptionType == ResDataTypes.Injection
    else
      return o.proSumptionType == ResDataTypes.Injection && !isSlack(o)
    end
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
        "sym_gen" => filter(x -> !isempty(x), [get_sym_gens(o) for o in net.prosumpsVec if filter_gens(o)]),
        "shunt" => filter(x -> !isempty(x), [get_shunts(o) for o in net.shuntVec]),
        "transformer" => filter(x -> !isempty(x), [get_trafos(t, net.baseMVA) for t in net.trafos if isa(t.comp, ImpPGMComp)]),
        "three_winding_transformer" => filter(x -> !isempty(x), [get_3WTtrafos(t, net.baseMVA) for t in net.trafos if isa(t.comp, ImpPGMComp3WT)]),
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
