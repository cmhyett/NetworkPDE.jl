using Test, NetworkPDE, ModelingToolkit, DomainSets

abstract type Pipe <: AbstractNetworkComponent end
abstract type FluxNode <: AbstractNetworkComponent end

function create_network()
    e1_params = Dict(:L => 100.0, :λ => 0.11, :D => 0.25, :from_node=>1, :to_node=>2, :type=>Pipe)
    e2_params = Dict(:L => 50.0, :λ => 0.11, :D => 0.25, :from_node=>2, :to_node=>3, :type=>Pipe)
    e1 = Edge(1, e1_params)
    e2 = Edge(2, e2_params)

    v1_params = Dict(:flux => zeros(3), :incident_edges=>[1], :type=>FluxNode);
    v2_params = Dict(:flux => zeros(3), :incident_edges=>[1,2], :type=>FluxNode);
    v3_params = Dict(:flux => zeros(3), :incident_edges=>[2], :type=>FluxNode);
    v1 = Vertex(1, v1_params)
    v2 = Vertex(2, v2_params)
    v3 = Vertex(3, v3_params)

    g = Network([e1, e2], [v1, v2, v3])
    return g;
end

@testset "constructing edge PDESystem" begin
    @parameters t
    net = create_network();
    @test get_metadata(net.edges[1], :type) == Pipe

    function instantiate(::Type{Pipe}, comp::AbstractNetworkComponent; name)
        @parameters x
        @variables ρ(..) ϕ(..)
        Dt = Differential(t)
        Dx = Differential(x)
        a = 1.0;
        eqs = [Dt(ρ(t,x)) + Dx(ϕ(t,x)) ~ 0,
              Dt(ϕ(t,x)) + a^2 * Dx(ρ(t,x)) ~ 0]
        ics = [ρ(0,x) ~ 0.0,#comp.initial_conditions[ρ]
               ϕ(0.0,x) ~ 0.0]#comp.initial_conditions[ϕ]]
        bcs = [] #fill in later
        domains = [x in Interval(0, get_metadata(comp, :L))]
        return PDESystem(eqs, [ics; bcs], domains, [t,x], [ρ(t,x), ϕ(t,x)], name=name);
    end
    @named sys = instantiate(get_metadata(net.edges[1], :type), net.edges[1]);
    @test sys isa PDESystem
end

@testset "compatibility with units" begin

end
