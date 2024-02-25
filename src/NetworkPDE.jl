module NetworkPDE

using Symbolics,
      OrdinaryDiffEq,
      ModelingToolkit,
      SciMLSensitivity,
      Zygote,
      Graphs,
      MetaGraphs,
      Plots,
      NetworkLayout

#@register_symbolic Base.floor(T::Type, x)::UInt64

include("graph_interface.jl")
#include("proof_of_concept.jl")

end
