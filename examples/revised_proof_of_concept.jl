using ModelingToolkit: t_nounits as t; #import global time var

function flux_node_factory!(v::Vertex)::PDESystem
    N = length(get_metadata(v, :incident_edges));
    function connect!(eqs, edges::Array{Edge, 1})
        #edge_sgns = 
    end
    sys = ODESystem(...);
    add_metadata(v, :sys, sys);
end

function edge_factory!(e::Edge, dxs)::PDESystem
    N = get_metadata(e, :L)/dxs[1]
    @parameters x
    @variables ρ(...), ϕ(...)
    eqs = [Dt(ρ(t,x)) + Dx(ϕ(t,x)) ~ 0,
           Dt(ϕ(t,x)) + a^2 * Dx(ρ(t,x)) ~ 0]
    bcs = [ρ(0, x) ~ e.initial_conditions[ρ],
           ϕ(0, x) ~ e.initial_conditions[ϕ]] # fill in later
    sys = PDESystem(...)
    add_metadata!(e, :sys, sys)
end

function create_network()::Network
    e1_params = Dict(:L => 100.0, :λ => 0.11, :D => 0.25, :from_node=1, :to_node=2, :instantiate=>edge_factory!)
    e2_params = Dict(:L => 50.0, :λ => 0.11, :D => 0.25, :from_node=>2, :to_node=>3, :instantiate=>edge_factory!)
    e1 = Edge(1, e1_params)
    e2 = Edge(2, e2_params)

    v1_params = Dict(:is_flux => true, :flux => zeros(3), :incident_edges=>[1], :instantiate=>flux_node_factory!)
    v2_params = Dict(:is_flux => true, :flux => zeros(3), :incident_edges=>[1,2], :instantiate=>flux_node_factory!)
    v3_params = Dict(:is_flux => true, :flux => zeros(3), :incident_edges=>[2], :instantiate=>flux_node_factory!)
    v1 = Vertex(1, v1_params)
    v2 = Vertex(2, v2_params)
    v3 = Vertex(3, v3_params)

    g = Network([e1, e2], [v1, v2, v3])
    return g;
end

network = create_network()
