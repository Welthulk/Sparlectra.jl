# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 03.05.2023
# include-file

queryLines = """
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX im: <http://imgpedia.dcc.uchile.cl/resource/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX entsoe: <http://entsoe.eu/CIM/SchemaExtension/3/1#>
PREFIX cim: <http://iec.ch/TC57/2013/CIM-schema-cim16#>
PREFIX md: <http://iec.ch/TC57/61970-552/ModelDescription/1#>

SELECT DISTINCT ?eqID ?lineName ?Unenn ?r ?r0 ?x ?x0 ?b ?b0 ?g ?g0 ?sCEndTemp ?length
WHERE {  
  ?Terminal a cim:Terminal ;                
            cim:Terminal.TopologicalNode ?topologicalNode ;
            cim:ACDCTerminal.connected "true" ;
            cim:Terminal.ConductingEquipment ?equipment .

  ?equipment a cim:ACLineSegment ;
             cim:ConductingEquipment.BaseVoltage/cim:BaseVoltage.nominalVoltage ?Unenn ;
             cim:IdentifiedObject.name ?lineName ;
             cim:Conductor.length ?length ;
             cim:ACLineSegment.r ?r ;
             cim:ACLineSegment.r0 ?r0 ;
             cim:ACLineSegment.x ?x ;
             cim:ACLineSegment.x0 ?x0 ;
             cim:ACLineSegment.bch ?b ;
             cim:ACLineSegment.b0ch ?b0 ;
             cim:ACLineSegment.gch ?g ;
             cim:ACLineSegment.g0ch ?g0 ;
             cim:ACLineSegment.shortCircuitEndTemperature ?sCEndTemp .

  BIND (strafter(str(?equipment), "#") AS ?eqID) 
}
ORDER BY DESC(xsd:float(?Unenn))
"""
function QueryLines!(url::String, ACLines::ACLineVector, debug::Bool)::Bool
  try
    response = HTTP.request("POST", url, headers = header2, body = queryLines)
    art = 1
    if HTTP.status(response) == 200

      result = JSON.parse(String(response.body))
      bindings = result["results"]["bindings"]
      if !isnothing(bindings)
        anz = length(bindings)
      else
        anz = 0
      end

      if anz == 0
        @warn "No lines found"
        return false
      end

      for binding in bindings
        ID = binding["eqID"]["value"]
        linename = binding["lineName"]["value"]
        Unenn = binding["Unenn"]["value"]
        r = binding["r"]["value"]
        #r0 = binding["r0"]["value"]          
        x = binding["x"]["value"]
        #x0 =  binding["x0"]["value"]                     
        g = binding["g"]["value"]
        #g0 =  binding["g0"]["value"]                     
        b = binding["b"]["value"]
        #b0 =  binding["b0"]["value"]
        length = binding["length"]["value"]
        #sCEndTemp = binding["sCEndTemp"]["value"]          

        #if debug
        # println("ID= ",ID,", L= ",linename,", Unenn =",Unenn,", r= ",r,", x= ",x,", g= ",g,", b= ",b,", length=",length)             
        #end           
        acseg = ResDataTypes.ACLineSegment(ID, linename, parse(Float64, Unenn), parse(Float64, length), parse(Float64, r), parse(Float64, x), parse(Float64, g), parse(Float64, b))
        push!(ACLines, acseg)
      end
      return true
    else
      if debug
        println("Error: $(HTTP.status(response)) $(HTTP.body(response))")
      end
      return false
    end
  catch err
    @error err
    return false
  end

end
