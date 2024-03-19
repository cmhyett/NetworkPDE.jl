# struct NetworkDiscretizer <: AbstractEquationSystemDiscretization
#     dxs
#     time
#     network
# end

# function instantiate(comp::AbstractNetworkComponent)
#     constructor(get_metadata(comp, :type), comp)
# end
# # examples
# struct FluxNode <: AbstractNetworkComponent
#     incoming_edges::Array{Edge, 1}
#     outgoing_edges::Array{Edge, 1}
#     params;
#     eqs;
# end

# struct Valve <: AbstractNetworkComponent
#     params;
#     eqs;
# end
