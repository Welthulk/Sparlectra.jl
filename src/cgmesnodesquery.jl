# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 28.04.2023
# include-file

queryNodes = """
PREFIX im: <http://imgpedia.dcc.uchile.cl/resource/>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX entsoe: <http://entsoe.eu/CIM/SchemaExtension/3/1#>
PREFIX md: <http://iec.ch/TC57/61970-552/ModelDescription/1#>
PREFIX cgmes: <http://iec.ch/TC57/2013/CIM-schema-cim16#>
PREFIX cim: <http://iec.ch/TC57/2013/CIM-schema-cim16#>

SELECT ?Vn ?nname ?NodeID ?tname ?eqname ?v ?angle ?p ?q ?shortClass ?eqID ?snum
WHERE {
  {
    ?Terminal a cim:Terminal ;
              cim:Terminal.TopologicalNode ?topologicalNode ;
              cim:Terminal.TopologicalNode/cim:TopologicalNode.BaseVoltage/cim:BaseVoltage.nominalVoltage ?Vn ;
              cim:ACDCTerminal.sequenceNumber ?snum ;
              cim:ACDCTerminal.connected "true" ;
              cim:Terminal.ConductingEquipment ?ConductingEquipment .

    OPTIONAL { ?Terminal cim:IdentifiedObject.name ?tname }

    ?ConductingEquipment rdf:type ?class ;
                         cim:IdentifiedObject.name ?eqname .

    ?topologicalNode cim:IdentifiedObject.name ?nodename ;
                     cim:IdentifiedObject.description ?desc .

    OPTIONAL {
      ?SvPowerFlow a cim:SvPowerFlow ;
                   cim:SvPowerFlow.Terminal ?Terminal ;
                   cim:SvPowerFlow.p ?p ;
                   cim:SvPowerFlow.q ?q .
    }

    OPTIONAL {
      ?SvVoltage a cim:SvVoltage ;
                 cim:SvVoltage.TopologicalNode ?topologicalNode ;
                 cim:SvVoltage.v ?v ;
                 cim:SvVoltage.angle ?angle .
    }

    FILTER (BOUND(?Terminal) && BOUND(?topologicalNode) && BOUND(?class))

    BIND (CONCAT(?nodename, "_", ?desc) AS ?nname)
    BIND (STRAFTER(STR(?ConductingEquipment), "#") AS ?eqID)
    BIND (STRAFTER(STR(?topologicalNode), "#") AS ?NodeID)
    BIND (STRAFTER(STR(?class), "#") AS ?shortClass)
  }
}
ORDER BY DESC(xsd:float(?Vn)) ASC(?NodeID) ASC(?eqname)
"""


function QueryNodes!(url::String, NodeVec::NodeVector, debug::Bool)::Bool
  try
    lastNodeName = ""
    lastU = ""
    lastNodeID = ""
    lastVm = nothing
    lastVa = nothing

    first_iteration = true
    terminals = ResDataTypes.Terminal[]

    function speichern()
      zone = nothing # TODO: Zone auslesen
      area = nothing # TODE: Area auslesen
      Vn = parse(Float64, lastU)

      if !isnothing(lastVm)
        vm_pu = lastVm / Vn
      else
        vm_pu = nothing
      end
      node = nothing
      try

        node = ResDataTypes.Node(lastNodeID, lastNodeName, Vn, copy(terminals))

      catch err
        @error err
        return false
      end
      if debug
        println("***")
        println("push:", node)
      end
      push!(NodeVec, node)
      empty!(terminals)
    end

    response = HTTP.request("POST", url, headers = header2, body = queryNodes)

    if HTTP.status(response) == 200

      result = JSON.parse(String(response.body))
      bindings = result["results"]["bindings"]
      if !isnothing(bindings)
        anz = length(bindings)
      else
        anz = 0
      end
      if anz == 0
        @warn "No nodes found"
        return false
      end
      # Loop over the query results
      #SELECT ?Vn ?nodename ?NodeID ?tname ?v ?angle ?p ?q  ?shortClass ?eqID ?snum
      for binding in bindings
        nodename = binding["nname"]["value"]
        Nennspannung = binding["Vn"]["value"]
        NodeID = binding["NodeID"]["value"]
        Vm = getFloatValue(binding, "v")
        Va = rad2deg(getFloatValue(binding, "angle"))

        #neuer Node, alles wegspeichern....
        if (lastNodeID != NodeID) || (first_iteration)
          #if (nodename != lastNodeName) || (first_iteration)
          if !first_iteration
            speichern()
          end
          first_iteration = false
          lastNodeName = nodename
          lastU = Nennspannung
          lastNodeID = NodeID
          lastVm = Vm
          lastVa = Va
        end

        pFlow = getFloatValue(binding, "p")
        qFlow = getFloatValue(binding, "q")

        TerminalName = binding["tname"]["value"]
        EquipmentName = binding["eqname"]["value"]

        Seite = binding["snum"]["value"]
        Art = binding["shortClass"]["value"]
        ID = binding["eqID"]["value"]
        if debug
          println("Vn= ", Nennspannung, ", node= ", nodename, ", Equipment= ", EquipmentName, " Art= ", Art, " ID = ", ID, " Seite= ", Seite, " p= ", pFlow, " q= ", qFlow)
        end

        cType = ResDataTypes.toComponentTyp(Art)

        if cType == ResDataTypes.UnknownC
          @warn "Unknown component type: ", Art
          continue
        end

        component = ResDataTypes.Component(ID, EquipmentName, ResDataTypes.toComponentTyp(Art), parse(Float64, Nennspannung))
        terminal = ResDataTypes.Terminal(component, ResDataTypes.toSeitenTyp(Seite))
        push!(terminals, terminal)

      end
      speichern()
      return true
    else
      if debug
        println("Error: $(HTTP.status(response)) $(HTTP.body(response))")
      end
      return false
    end
  catch err

    #st= stacktrace()
    #bt = catch_backtrace()[1]
    @error err

    return false
  end
end
