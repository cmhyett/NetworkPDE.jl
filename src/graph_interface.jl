abstract type AbstractNetworkComponent end

mutable struct Edge <: AbstractNetworkComponent
    id::Int
    metadata::Any
end
Edge(id) = Edge(id, Dict());

mutable struct Vertex <: AbstractNetworkComponent
    id::Int
    metadata::Any
end
Vertex(id) = Vertex(id, Dict());

mutable struct Network <: AbstractNetworkComponent
    edges::Array{Edge}
    vertices::Array{Vertex}
    metadata::Any
end
Network() = Network(Edge[], Vertex[], Dict());
Network(edges, vertices) = Network(edges, vertices, Dict());

add_metadata(a::AbstractNetworkComponent, key, data) = a.metadata[key] = data;
function add_metadata(a::AbstractNetworkComponent, d::Dict...)
    (a.metadata = merge(a.metadata, d...))
end;
get_metadata(a::AbstractNetworkComponent, key) = a.metadata[key];
add_edge(n::Network, e::Edge) = push!(n.edges, e);
get_edge(n::Network, id::Int) = n.edges[id];
get_vertex(n::Network, id::Int) = n.vertices[id];
num_edges(n::Network) = length(n.edges);
num_vertices(n::Network) = length(n.vertices);
