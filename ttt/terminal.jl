# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file termina.jl

# Data type to describe a terminal 
struct Terminal
  comp::AbstractComponent  
  
  function Terminal(comp::AbstractComponent)    
    new(comp)
  end

  function Base.show(io::IO, terminal::Terminal)
    print(io, "Terminal(")
    print(io, terminal.comp)    
  end
end
