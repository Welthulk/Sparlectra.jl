# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file termina.jl

# Data type to describe a terminal 
mutable struct Terminal
  comp::AbstractComponent  
  seite::SeitenTyp # Seite1 = von, Seite2 = zu, Seite3 = 3WT, tertiäre Seite

  function Terminal(comp::AbstractComponent, seite::SeitenTyp)
    new(comp, seite)
  end

  function Base.show(io::IO, terminal::Terminal)
    print(io, "Terminal(")
    print(io, terminal.comp)
    print(io, " Seite: ", terminal.seite, ")")
  end
end

function ChageSequnceNumber!(t::Terminal, n::SeitenTyp)
  t.seite = n
end
