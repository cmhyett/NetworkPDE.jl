using Test, NetworkPDE, ModelingToolkit, DomainSets

abstract type Pipe <: AbstractNetworkComponent end
abstract type FluxNode <: AbstractNetworkComponent end

function create_network()
    dx = 50
    e1_params = Dict(
        :L => 100.0, :λ => 0.0, :D => 0.25, :from_node => 1, :to_node => 2, :type => Pipe, :dx => dx)
    e2_params = Dict(
        :L => 150.0, :λ => 0.0, :D => 0.25, :from_node => 2, :to_node => 3, :type => Pipe, :dx => dx)
    e1 = Edge(1, e1_params)
    e2 = Edge(2, e2_params)

    v1_params = Dict(:flux => zeros(3), :incident_edges => [1], :type => FluxNode)
    v2_params = Dict(:flux => zeros(3), :incident_edges => [1, 2], :type => FluxNode)
    v3_params = Dict(:flux => zeros(3), :incident_edges => [2], :type => FluxNode)
    v1 = Vertex(1, v1_params)
    v2 = Vertex(2, v2_params)
    v3 = Vertex(3, v3_params)

    g = Network([e1, e2], [v1, v2, v3])
    return g
end

@testset "constructing edge PDESystem" begin
    @parameters t
    net = create_network()
    @test get_metadata(net.edges[1], :type) == Pipe

    function instantiate(::Type{Pipe}, comp::AbstractNetworkComponent; name)
        @parameters x
        @variables ρ(..) ϕ(..)
        Dt = Differential(t)
        Dx = Differential(x)
        a = 1.0
        eqs = [Dt(ρ(t, x)) + Dx(ϕ(t, x)) ~ 0,
            Dt(ϕ(t, x)) + a^2 * Dx(ρ(t, x)) ~ 0]
        ics = [ρ(0, x) ~ 0.0,#comp.initial_conditions[ρ]
            ϕ(0.0, x) ~ 0.0]#comp.initial_conditions[ϕ]]
        bcs = [] #fill in later
        domains = [x in Interval(0, get_metadata(comp, :L))]
        return PDESystem(eqs, [ics; bcs], domains, [t, x], [ρ(t, x), ϕ(t, x)], name = name)
    end
    @named sys = instantiate(get_metadata(net.edges[1], :type), net.edges[1])
    @test sys isa PDESystem
end

@testset "constructing node PDESystem" begin
    @parameters t
    net = create_network()
    @test get_metadata(net.vertices[1], :type) == FluxNode

    function instantiate(::Type{FluxNode}, comp::AbstractNetworkComponent; name)
        @parameters q, dx, dt
        @variables ρ(t)
        Dt = Differential(t)
        inc_edges = [net.edges[i] for i in get_metadata(comp, :incident_edges)]
        num_inc_edges = length(inc_edges)
        @variables tmp_ρ(t)[1:num_inc_edges] tmp_ϕ(t)[1:num_inc_edges]
        sgns = [get_metadata(e, :to_node) == comp.id ? 1 : -1 for e in inc_edges]
        cross_sections = [π * (get_metadata(e, :D) / 2)^2 for e in inc_edges]
        eqs = [ρ * (dx / dt) * sum(cross_sections) ~ q +
                                                     sum([(dx / dt) * cross_sections[i] *
                                                          tmp_ρ[i] -
                                                          sgns[i] * cross_sections[i] *
                                                          tmp_ϕ[i] for i in 1:num_inc_edges])]
        ics = []#ρ(0) ~ 0.0]#comp.initial_conditions[ρ]]
        bcs = [] #fill in later
        return ODESystem(eqs, t, name = name)
    end
    @named sys = instantiate(get_metadata(net.vertices[1], :type), net.vertices[1])
    @test sys isa ODESystem
end

@testset "completing dynamic network" begin
    Wave_Speed = 350
    @parameters t
    net = create_network()
    pressure_from_density(density) = Wave_Speed^2 * density
    density_from_pressure(pressure) = pressure / Wave_Speed^2
    #CFL ~ min(dx)/dt < sqrt(P'(ρ))
    dt = minimum([get_metadata(net.edges[i], :dx) for i in 1:length(net.edges)]) /
         sqrt(gradient(pressure_from_density, 1.0)[1])

    #    @test get_metadata(net.edges[1], :type) == Pipe

    function instantiate(::Type{Pipe}, comp::AbstractNetworkComponent; name)
        #        @parameters x
        Dt = Differential(t)
        #        Dx = Differential(x)

        L = get_metadata(comp, :L)
        dx = get_metadata(comp, :dx)
        β = get_metadata(comp, :λ) / (2 * get_metadata(comp, :D))
        nx = floor(Int, L / dx)
        @variables ρ(t)[1:nx] ϕ(t)[1:(nx + 1)] #staggered grid

        function dρ(ρ, ϕ, dx, dt)
            return @. -(1 / dx) * (ϕ[2:end] - ϕ[1:(end - 1)])
        end

        Finv(y, b) = y #sign.(y) .* ((-1 .+ sqrt.(1 .+ (4*b).*abs.(y)))./(2b))

        function dϕ(ρ, ϕ, dx, dt, β, pressure_from_density)
            b = (β * dt) ./ (ρ[1:(end - 1)] + ρ[2:end])
            y = ϕ[2:(end - 1)] .-
                (dt / dx) .* (pressure_from_density(ρ[2:end]) .-
                 pressure_from_density(ρ[1:(end - 1)])) .-
                b .* (ϕ[2:(end - 1)] .* abs.(ϕ[2:(end - 1)]))
            return (Finv(y, b) .- ϕ[2:(end - 1)]) ./ dt
        end

        drho = dρ(ρ, ϕ, dx, dt)
        internal_eqs = vcat(Symbolics.scalarize.([Dt(ρ) ~ drho;
                                                  Dt(ϕ[2:(end - 1)]) ~ dϕ(ρ, ϕ, dx, dt, β,
                                                      pressure_from_density)])...)
        ics = [ρ .~ 0.0,#comp.initial_conditions[ρ]
            ϕ .~ 0.0]#comp.initial_conditions[ϕ]]
        bcs = [] #fill in later
        #domains = [x in Interval(0, L)]

        #return PDESystem(eqs, [ics; bcs], domains, [t, x], [ρ(t, x), ϕ(t, x)], name = name)
        return ODESystem(internal_eqs, t, name = name)
    end

    function instantiate(::Type{FluxNode}, comp::AbstractNetworkComponent; name)
        @parameters q
        @variables ρ(t)

        return ODESystem([0 ~ 0], t, [ρ], [q], name = name)
    end

    function connect(net::Network; name)
        coupling_eqs = []
        for v in net.vertices
            inc_edges = [net.edges[i] for i in get_metadata(v, :incident_edges)]
            num_inc_edges = length(inc_edges)
            sgns = [get_metadata(e, :to_node) == v.id ? 1 : -1 for e in inc_edges]
            cross_sections = [π * (get_metadata(e, :D) / 2)^2 for e in inc_edges]
            edge_dxs = [get_metadata(e, :dx) for e in inc_edges]
            edge_ρs = []
            edge_ϕs = []
            node_ρ = unknowns(v.metadata[:model], v.metadata[:model].var_to_name[:ρ])
            node_q = unknowns(v.metadata[:model], v.metadata[:model].var_to_name[:q])
            for i in 1:num_inc_edges
                model = inc_edges[i].metadata[:model]
                if (sgns[i] == 1)
                    push!(edge_ρs,
                        unknowns(
                            model, model.var_to_name[:ρ][length(model.var_to_name[:ρ])]))
                    push!(edge_ϕs,
                        unknowns(
                            model, model.var_to_name[:ϕ][length(model.var_to_name[:ϕ])]))
                else
                    push!(edge_ρs, unknowns(model, model.var_to_name[:ρ][1]))
                    push!(edge_ϕs, unknowns(model, model.var_to_name[:ϕ][1]))
                end
            end
            push!(coupling_eqs,
                node_ρ * (1 / dt) * sum([edge_dxs[i] * cross_sections[i] for i in 1:num_inc_edges]) ~ node_q +
                                                                                                      sum([(edge_dxs[i] /
                                                                                                            dt) *
                                                                                                           cross_sections[i] *
                                                                                                           edge_ρs[i] -
                                                                                                           sgns[i] *
                                                                                                           cross_sections[i] *
                                                                                                           edge_ϕs[i]
                                                                                                           for i in 1:num_inc_edges]))
            for i in 1:num_inc_edges
                push!(coupling_eqs,
                    Dt(edge_ϕs[i]) ~ sgns[i] * (pressure_from_density(node_ρ) -
                                      pressure_from_density(edge_ρs[i])))
            end
        end
        return compose(ODESystem(coupling_eqs, t, name = name),
            [get_metadata(net.edges[i], :model) for i in 1:length(net.edges)]...,
            [get_metadata(net.vertices[i], :model) for i in 1:length(net.vertices)]...)
    end

    function instantiate(network::Network; name)
        connection_eqs = []
        for i in 1:length(network.edges)
            add_metadata(network.edges[i], :model,
                instantiate(get_metadata(network.edges[i], :type),
                    network.edges[i], name = Symbol("edge_$(i)")))
        end
        for i in 1:length(network.vertices)
            add_metadata(network.vertices[i], :model,
                instantiate(get_metadata(network.vertices[i], :type),
                    network.vertices[i], name = Symbol("node_$(i)")))
        end
        return connect(network, name = name)
    end
    @named sys = instantiate(net)
    @test sys isa ODESystem
end

@testset "compatibility with units" begin end
