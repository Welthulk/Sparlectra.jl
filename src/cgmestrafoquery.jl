# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 19.05.2023
# include-file
#

queryTrafos = """
PREFIX im: <http://imgpedia.dcc.uchile.cl/resource/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX entsoe: <http://entsoe.eu/CIM/SchemaExtension/3/1#>
PREFIX cim: <http://iec.ch/TC57/2013/CIM-schema-cim16#>
PREFIX md: <http://iec.ch/TC57/61970-552/ModelDescription/1#>

SELECT DISTINCT ?endNumber ?name ?eqID ?Unenn ?r ?x ?b ?g ?ratedU ?ratedS ?GeneratorTrafo ?NameStufenSteller ?enable ?Stellung ?Unten ?Oben ?Mittelstellung ?DefaultStellung ?stepV ?nU
WHERE {
  ?Terminal a cim:Terminal ;
            cim:Terminal.TopologicalNode ?topologicalNode ;
            cim:ACDCTerminal.connected "true" ;
            cim:Terminal.ConductingEquipment ?equipment .

  ?topologicalNode cim:IdentifiedObject.name ?nodename ;
                   cim:IdentifiedObject.description ?desc .
  ?equipment a cim:PowerTransformer ;
             cim:IdentifiedObject.name ?name .

  ?TransformerEnd a cim:PowerTransformerEnd ;
                 cim:PowerTransformerEnd.PowerTransformer ?equipment ;
                 cim:TransformerEnd.BaseVoltage/cim:BaseVoltage.nominalVoltage ?Unenn ;
                 cim:PowerTransformerEnd.b ?b ;
                 cim:PowerTransformerEnd.r ?r ;
                 cim:PowerTransformerEnd.x ?x ;
                 cim:PowerTransformerEnd.ratedU ?ratedU ;
                 cim:PowerTransformerEnd.ratedS ?ratedS ;
                 cim:TransformerEnd.endNumber ?endNumber ;
                 cim:TransformerEnd.Terminal ?terminal .

  OPTIONAL {
    ?equipment cim:Equipment.EquipmentContainer ?EquipmentContainer .
    ?equipment cim:PowerTransformer.isPartOfGeneratorUnit ?GeneratorTrafo .
    ?TransformerEnd cim:PowerTransformerEnd.connectionKind ?connectionKind .
    ?TransformerEnd cim:PowerTransformerEnd.g ?g .
    ?RatioTapChanger a cim:RatioTapChanger ;
                     cim:IdentifiedObject.name ?NameStufenSteller ;
                     cim:TapChanger.controlEnabled ?enable ;
                     cim:TapChanger.step ?Stellung ;
                     cim:TapChanger.lowStep ?Unten ;
                     cim:TapChanger.highStep ?Oben ;
                     cim:TapChanger.neutralStep ?Mittelstellung ;
                     cim:TapChanger.normalStep ?DefaultStellung ;
                     cim:TapChanger.neutralU ?nU ;
                     cim:RatioTapChanger.stepVoltageIncrement ?stepV ;
                     cim:RatioTapChanger.TransformerEnd ?TransformerEnd .
  }

  BIND (strafter(str(?equipment), "#") AS ?eqID) .
}
ORDER BY ASC(?name) ASC(?endNumber)
"""

function QueryTrafos!(url::String, trafos::TrafoVector, debug::Bool)::Bool
  try

    response = HTTP.request("POST", url, headers = header2, body = queryTrafos)

    if HTTP.status(response) == 200

      result = JSON.parse(String(response.body))
      bindings = result["results"]["bindings"]
      anz = 0
      if !isnothing(bindings)
        anz = length(bindings)
      else
        anz = 0
      end
      if anz == 0
        @warn "No transformers found"
        return false
      end

      len = length(result)

      # erstes Element
      t = bindings[1]
      lastName = t["name"]["value"]

      rangeVec = []
      change = false
      r1 = 0
      r2 = 0
      firstIndex = 1
      lastIndex = 0
      # baut einen Vector mit Bereichen gleicher Trafoseiten auf 1:2, 3:5, .....     
      for i in eachindex(bindings)
        t = bindings[i]
        N = t["name"]["value"]

        if change
          lastName = N

          firstIndex = i
          change = false
        end

        if i == length(bindings)
          r1 = r2 + 1
          r2 = length(bindings)

          push!(rangeVec, (r1:r2))
        elseif N != lastName && i != length(bindings)
          if firstIndex == 1
            r1 = firstIndex
          else
            r1 = firstIndex - 1
          end
          r2 = lastIndex
          change = true

          push!(rangeVec, (r1:r2))
        else
          lastIndex = i

        end
      end

      if debug
        for i in rangeVec
          @show i
        end
      end

      i = 0
      shift_degree = 0.0
      # durchlaufe den Vector mit den jeweiligen Bereichen
      for k in rangeVec
        cmp = nothing
        s1 = nothing
        s2 = nothing
        s3 = nothing
        tap = nothing
        bereich = k
        anzTaps = 0
        len = length(bereich)
        for j in bereich
          i += 1
          regelungEin = false
          binding = bindings[j]
          ID = binding["eqID"]["value"]
          Name = binding["name"]["value"]
          Unenn = getFloatValue(binding, "Unenn")
          Seite = getIntegerValue(binding, "endNumber")
          ratedU = getFloatValue(binding, "ratedU")
          ratedS = getFloatValue(binding, "ratedS")
          r = getFloatValue(binding, "r")
          x = getFloatValue(binding, "x")
          b = getFloatValue(binding, "b")
          g = getFloatValue(binding, "g")

          if debug
            @show ID, Name, Unenn, Seite, ratedU, ratedS, r, x, b, g
          end

          nameStufenSteller = getStringValue(binding, "NameStufenSteller")

          if !(isStringEmpty(nameStufenSteller))
            if debug
              println("StufenSteller: ", nameStufenSteller)
            end
            anzTaps += 1
            u = getIntegerValue(binding, "Unten")
            o = getIntegerValue(binding, "Oben")
            m = getIntegerValue(binding, "Mittelstellung")
            n = getIntegerValue(binding, "DefaultStellung")
            stl = getIntegerValue(binding, "Stellung")
            regelungEin = getBoolValue(binding, "enable")


            nV = getFloatValue(binding, "nU")
            incV = getFloatValue(binding, "stepV")

            if debug
              println("tap $(nameStufenSteller)")
              @show u, o, m, n, stl, regelungEin, incV, nV
              @show stl, u, o, m, n, nV, incV
            end
            tap = PowerTransformerTaps(stl, u, o, m, incV, nV)
            if debug
              @show tap
            end
          else
            tap = nothing
          end

          if Seite == 1
            cmp = ResDataTypes.Component(ID, Name, "POWERTRANSFORMER", Unenn)
            if debug
              @show cmp
            end
            s1 = ResDataTypes.PowerTransformerWinding(Unenn, r, x, b, g, shift_degree, ratedU, ratedS, tap)
            if debug
              @show "Seite 1 =", s1
            end
          elseif Seite == 2
            s2 = ResDataTypes.PowerTransformerWinding(Unenn, r, x, b, g, shift_degree, ratedU, ratedS, tap)

            if debug
              @show "Seite 2 =", s2
            end

          elseif Seite == 3
            s3 = ResDataTypes.PowerTransformerWinding(Unenn, r, x, b, g, shift_degree, ratedU, ratedS, tap)

            if debug
              @show "Seite 3 =", s3
            end
          end
          # Ende des Bereichs erreicht, abspeichern....    
          if j == bereich[len]
            Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
            push!(trafos, Trafo)
            if debug
              @show "Trafo:", Trafo
            end
          end
        end
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

end # QueryTrafos!


