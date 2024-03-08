# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 08.06.2023
# include-file cgmesprosumptionquery.jl
# purpose: Query for producers and consumptions data
#

#=
Example ResultSet:
"Unenn","NodeName","NodeID","eqName","Snummer","nodeType","NodeType","shortClass","eqID","ratedS","ratedU","maxP","maxQ","minP","minQ","referencePriority","ratedPowerFactor"
"380","1","_adee76cd-b2b9-48ac-8fd4-0d205a435f59","Q1","1","PQ","E","ExternalNetworkInjection","_089c1945-4101-487f-a557-66c013b748f6","","","800","600","-800","-600","0",""
"110","5","_37edd845-456f-4c3e-98d5-19af0c1cef1e","Q2","1","PQ","E","ExternalNetworkInjection","_3de9e1ad-4562-44df-b268-70ed0517e9e7","","","88","66","-88","-66","0",""
"21","HG1","_7f5515b2-ca6b-45af-93ee-f196686f0c66","G1","1","PV","E","SynchronousMachine","_ca67be42-750e-4ebf-bfaa-24d446e59a22","150","21","","79","","-79","0","0.85"
"10","6","_764e0b8a-f2af-4092-b6aa-b4a19e55db98","G3","1","PV","E","SynchronousMachine","_392ea173-4f8e-48fa-b2a3-5c3721e93196","10","10.5","","6","","-6","0","0.8"
"10","7","_cd84fa40-ef63-422d-8ee0-d0a0f806719e","M2a","1","PQ","L","AsynchronousMachine","_ba62884d-8800-41a8-9c26-698297d7ebaa","2.321","10","","","","","","0.89"
"10","7","_cd84fa40-ef63-422d-8ee0-d0a0f806719e","M2b","1","PQ","L","AsynchronousMachine","_f184d87b-5565-45ee-89b4-29e8a42d3ad1","2.321","10","","","","","","0.89"
"10","7","_cd84fa40-ef63-422d-8ee0-d0a0f806719e","M3","1","PQ","L","AsynchronousMachine","_062ece1f-ade5-4d20-9c3a-fd8f12d12ec1","5.828","10","","","","","","0.88"
"10","HG2","_c7eda3d2-e92d-4935-8166-5e045d3de045","G2","1","PV","E","SynchronousMachine","_2970a2b7-b840-4e9c-b405-0cb854cd2318","100","10.5","","43.6","","-43.6","1","0.9"
=#



QueryProSumptionStr = """
PREFIX im: <http://imgpedia.dcc.uchile.cl/resource/>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX entsoe: <http://entsoe.eu/CIM/SchemaExtension/3/1#>
PREFIX cim: <http://iec.ch/TC57/2013/CIM-schema-cim16#>
PREFIX md: <http://iec.ch/TC57/61970-552/ModelDescription/1#>
PREFIX math: <http://www.w3.org/2005/xpath-functions/math#>

SELECT ?Vn ?NodeID ?eqID ?eqName ?shortClass ?pVal ?qVal ?minP ?maxP ?minQ ?maxQ ?qPercent ?ratedS ?ratedPF ?ratedU ?referencePriority
WHERE 
{       
    ?Terminal a cim:Terminal ;                
              cim:Terminal.TopologicalNode ?topologicalNode;
              cim:Terminal.TopologicalNode/cim:TopologicalNode.BaseVoltage/cim:BaseVoltage.nominalVoltage ?Vn;
              cim:ACDCTerminal.sequenceNumber ?snum;
              cim:ACDCTerminal.connected "true";                
              cim:Terminal.ConductingEquipment ?EnergyConsumer.
    ?EnergyConsumer rdf:type ?type;
                    cim:IdentifiedObject.name ?eqName.   
    VALUES ?type { cim:SynchronousMachine cim:EnergyConsumer cim:ExternalNetworkInjection cim:AsynchronousMachine cim:GeneratingUnit }
    
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.qPercent ?qPercent }
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.minP ?minPval }    
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.minQ ?minQval }    
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.maxP ?maxPval }    
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.maxQ ?maxQval }
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.p  ?pEVal}
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.q  ?qEVal}
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.minP ?minEP}
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.maxP ?maxEP}
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.minQ ?minEQ}
    OPTIONAL { ?EnergyConsumer cim:ExternalNetworkInjection.maxQ ?maxEQ}
  
    OPTIONAL { ?EnergyConsumer cim:RotatingMachine.p ?pValue}
    OPTIONAL { ?EnergyConsumer cim:RotatingMachine.q ?qValue}
  
    OPTIONAL { ?EnergyConsumer cim:RotatingMachine.ratedS ?ratedS }    
    OPTIONAL { ?EnergyConsumer cim:RotatingMachine.ratedU ?ratedU }
    OPTIONAL { ?EnergyConsumer cim:RotatingMachine.ratedPowerFactor ?ratedPF }
    OPTIONAL { ?EnergyConsumer cim:SynchronousMachine.referencePriority ?referencePriority }
    
    BIND (strafter(str(?EnergyConsumer), "#") AS ?eqID)  
    BIND (strafter(str(?type), "#") AS ?shortClass)
    BIND (strafter(str(?topologicalNode), "#") AS ?NodeID)
    
    BIND (COALESCE(?pValue, ?pEVal) AS ?pVal)
    BIND (COALESCE(?qValue, ?qEVal) AS ?qVal)
    BIND (COALESCE(?minPval, ?minEP) AS ?minP)
    BIND (COALESCE(?maxPval, ?maxEP) AS ?maxP)
    BIND (COALESCE(?minQval, ?minEQ) AS ?minQ)    
    BIND (COALESCE(?maxQval, ?maxEQ) AS ?maxQ)  
}
ORDER BY DESC(xsd:float(?Vn))
"""

function QueryProSumption!(url::String, prosumptionVec::ProSumptionVector, debug::Bool)::Bool
  try

    response = HTTP.request("POST", url, headers = header2, body = QueryProSumptionStr)

    if HTTP.status(response) == 200

      result = JSON.parse(String(response.body))
      bindings = result["results"]["bindings"]
      if !isnothing(bindings)
        anz = length(bindings)
      else
        anz = 0
      end
      if anz == 0
        @warn "No prosumers found"
        return false
      end

      #SELECT ?Vn ?NodeID ?eqID ?shortClass ?p ?q ?minP ?minQ ?maxP ?maxQ ?qPercent ?ratedS ?ratedPF ?ratedU ?referencePriority
      for binding in bindings
        shortClass = binding["shortClass"]["value"]
        Vn = getFloatValue(binding, "Vn")
        eqID = binding["eqID"]["value"]
        nodeID = binding["NodeID"]["value"]
        eqName = binding["eqName"]["value"]

        qPercent = getFloatValue(binding, "qPercent")
        pVal = getFloatValue(binding, "pVal")
        qVal = getFloatValue(binding, "qVal")
        minP = getFloatValue(binding, "minP")
        minQ = getFloatValue(binding, "minQ")
        maxP = getFloatValue(binding, "maxP")
        maxQ = getFloatValue(binding, "maxQ")

        ratedS = getFloatValue(binding, "ratedS")
        ratedU = getFloatValue(binding, "ratedU")
        ratedPF = getFloatValue(binding, "ratedPF")
        referencePriority = getIntegerValue(binding, "referencePriority")

        sClass = ResDataTypes.toComponentTyp(shortClass)
        c = ResDataTypes.Component(eqID, eqName, sClass, Vn)

        o = ResDataTypes.ProSumer(c, nodeID, ratedS, ratedU, qPercent, pVal, qVal, maxP, minP, maxQ, minQ, ratedPF, referencePriority)
        if debug
          @show o
        end
        push!(prosumptionVec, o)
      end
      return true
    else
      println("Error: ", HTTP.status(response))
      return false
    end
  catch e
    @error e
    return false
  end

end