abstract type AbstractNetworkComponent end
    
mutable struct Edge <: AbstractNetworkComponent
    id::Int;
    metadata;
end
Edge(id) = Edge(id, Dict());

mutable struct Vertex <: AbstractNetworkComponent
    id::Int;
    metadata;
end
Vertex(id) = Vertex(id, Dict());

mutable struct Network <: AbstractNetworkComponent
    edges::Array{Edge};
    vertices::Array{Vertex};
end
Network() = Network(Edge[], Vertex[]);

add_metadata(a::AbstractNetworkComponent, key, data) = a.metadata[key] = data;
add_metadata(a::AbstractNetworkComponent, d::Dict...) = merge(a.metadata, d...);
get_metadata(a::AbstractNetworkComponent, key) = a.metadata[key];
add_edge(
get_edge(n::Network, id::Int) = n.edges[id];
get_vertex(n::Network, id::Int) = n.vertices[id];
