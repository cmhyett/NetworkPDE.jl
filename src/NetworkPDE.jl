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
include("problem_interface.jl")
#include("proof_of_concept.jl")
export AbstractNetworkComponent, Edge, Vertex, Network
export add_metadata, get_metadata, add_edge, get_edge, get_vertex, num_edges, num_vertices

export AbstractNetworkComponent, Edge, Vertex, Network, add_metadata, get_metadata, add_edge, get_edge, get_vertex, num_edges, num_vertices

end
