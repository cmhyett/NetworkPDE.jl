
function graphplot(g::AbstractGraph)
    p = plot()
    v_positions = sfdp(g)
    for e in edges(g)
        plot!(p,
              [v_positions[e.src][1], v_positions[e.dst][1]],
              [v_positions[e.src][2], v_positions[e.dst][2]],
              color = :black, linewidth=5)
    end
    scatter!(p, v_positions, markersize=8, color=:blue)
    plot!(p, legend=false)
    return p
end


# interior updates
function dρ(ρ, ϕ, dx, dt)
    return @. -(1/dx)*(ϕ[2:end]-ϕ[1:end-1])
end
function dϕ_interior(ρ, ϕ, dx, dt)
    return @. -(dt/dx) * (p(ρ[2:end]) - p(ρ[1:end-1])) - (β*dt)*(ϕ[2:end-1]*abs(ϕ[2:end-1]))/(ρ[1:end-1]+ρ[2:end])
end
function get_dvs(g::MetaGraph)
    edge_dvs = vcat([[get_prop(g, edge, :dvs)...] for edge in edges(g)]...);
    vertex_dvs = [get_prop(g, vertex, :dvs) for vertex in vertices(g)];
    return vcat(vertex_dvs..., edge_dvs...);
end
function get_params(g::MetaGraph)
    params = [];
    for edge in edges(g)
        has_prop(g, edge, :params) ? push!(params, [get_prop(g, edge, :params)...]) : continue;
    end
    for vertex in vertices(g)
        has_prop(g, vertex, :params) ? push!(params, [get_prop(g, vertex, :params)...]) : continue;
    end
    return vcat(params...);
end

function create_system(g::AbstractGraph, t, dx, params)
    i = 0;
    eqs = [];
    dx = 25.0;
    Dt = Differential(t)
    
    
    function get_incoming_edges(vertex, graph::MetaGraph)::Vector{Graphs.SimpleGraphs.SimpleEdge}
        edge_list = [e.dst == vertex ? e : nothing for e in edges(g)];
        return edge_list[findall(x->x!=nothing, edge_list)]
    end
    function get_outgoing_edges(vertex, graph::MetaGraph)::Vector{Graphs.SimpleGraphs.SimpleEdge}
        edge_list = [e.src == vertex ? e : nothing for e in edges(g)];
        return edge_list[findall(x->x!=nothing, edge_list)]
    end
    function get_adj_edges(vertex, graph::MetaGraph)::Vector{Graphs.SimpleGraphs.SimpleEdge}
        edge_list = [(e.src == vertex || e.dst == vertex) ? e : nothing for e in edges(g)];
        return edge_list[findall(x->x!=nothing, edge_list)]
    end
    function get_edge_signs(vertex, edges::Vector{Graphs.SimpleGraphs.SimpleEdge})::Vector{Int}
        return [e.dst == vertex ? 1 : -1 for e in edges];
    end
    function get_cross_section(edge)
        return π*(get_prop(g, edge, :D)/2)^2;
    end
    function get_bdy_flux(vertex, params, t)
        return params[(vertex-1)*length(params_tsteps)+floor(Int, t/param_dt)+1];
    end


    for edge in edges(g)
        i += 1;
        xs = 0:dx:get_prop(g, edge, :L)
        num_dx = length(xs);
        sym = Symbol(":ρ_e$(i)");
        ρ = only(@variables($sym(t)[1:num_dx]))
        sym = Symbol(":ϕ_e$(i)");
        ϕ = only(@variables($sym(t)[1:num_dx+1]))
        append!(eqs, Symbolics.scalarize.([Dt.(ρ) .~ dρ(ρ, ϕ, dx, dt);
                                           Dt.(ϕ[2:end-1]) .~ dϕ_interior(ρ, ϕ, dx, dt)]))
        set_prop!(g, edge, :dvs, [ρ, ϕ])
    end

    for vertex in vertices(g)
        adj_edges = get_adj_edges(vertex, g);
        edge_sgns = Dict(adj_edges .=> get_edge_signs(vertex, adj_edges));
        sym = Symbol(":ρ_n$(vertex)");
        ρ = only(@variables($sym(t)))
        sym = Symbol(":ϕ_n$(vertex)");
        ϕ = only(@variables($sym(t)))
        sym = Symbol(":p_n$(vertex)");
        #params = only(@parameters($sym[1:length(params_tsteps)]))
        q = get_bdy_flux(vertex, params, t)
        rhs_terms = [q]
        lhs_terms = []
        for edge in adj_edges
            e_ρ, e_ϕ = get_prop(g, edge, :dvs);
            cross_section = get_cross_section(edge);
            if (edge_sgns[edge] == 1) #outgoing
                push!(rhs_terms, (dx/dt)*cross_section*e_ρ[1]);
                push!(rhs_terms, -edge_sgns[edge]*cross_section*e_ϕ[1])
            else #incoming
                push!(rhs_terms, (dx/dt)*cross_section*e_ρ[end]);
                push!(rhs_terms, edge_sgns[edge]*cross_section*e_ϕ[end])
            end
            push!(lhs_terms, (dx/dt)*cross_section)
        end
        append!(eqs, [Dt(ρ) ~ ((sum(rhs_terms)/sum(lhs_terms) - ρ)/dt)])
        set_prop!(g, vertex, :dvs, [ρ, ϕ])
#        set_prop!(g, vertex, :params, params)
    end

    function dϕ_bdy(ρ, ϕ, dx, dt) #length(ρ)==2,length(ϕ)=1
        return -(dt/dx) * (p(ρ[2]) - p(ρ[1])) - (β*dt)*(ϕ[1]*abs(ϕ[1]))/(ρ[1]+ρ[2])
    end
    function dϕ_lbdy(g, edge, dx, dt)
        edge_ρ, edge_ϕ = get_prop(g, edge, :dvs);
        node_ρ, node_ϕ = get_prop(g, edge.src, :dvs);
        return dϕ_bdy([node_ρ; edge_ρ[1]], edge_ϕ[1], dx, dt)
    end
    function dϕ_rbdy(g, edge, dx, dt)
        edge_ρ, edge_ϕ = get_prop(g, edge, :dvs);
        node_ρ, node_ϕ = get_prop(g, edge.src, :dvs);
        return dϕ_bdy([edge_ρ[end]; node_ρ], edge_ϕ[end], dx, dt)
    end
    # finally, now that vertices have variables defined, finish interior discretization
    for edge in edges(g)
        ρ, ϕ = get_prop(g, edge, :dvs);
        append!(eqs, [Dt(ϕ[1]) ~ dϕ_lbdy(g, edge, dx, dt);
                      Dt(ϕ[end]) ~ dϕ_rbdy(g, edge, dx, dt)]);
    end
    eqs = Vector{Equation}(eqs);
    dvs = get_dvs(g);
#    params = [(get_prop(g, v, :params)).value for v in vertices(g)]#get_params(g);
    @named de = ODESystem(eqs, t, dvs, [params...])
    return de
end

function run_example()
    # specify network
    g = MetaGraph()
    v1_params = Dict(:is_flux => true, :flux => zeros(3))
    v2_params = Dict(:is_flux => true, :flux => zeros(3))
    v3_params = Dict(:is_flux => true, :flux => zeros(3))
    vertex_dict = Dict(1 => v1_params, 2 => v2_params, 3 => v3_params)
    e1_params = Dict(:L => 100.0, :λ => 0.11, :D => 0.25)
    e2_params = Dict(:L => 50.0, :λ => 0.11, :D => 0.25)
    edge_dict = Dict(1 => e1_params, 2 => e2_params)
    [add_vertex!(g, vertex_dict[i]) for i = 1:length(vertex_dict)]
    add_edge!(g, 1, 2, e1_params)
    add_edge!(g, 2, 3, e2_params)

    # discretize network
    @variables t
    @parameters params[1:9] #current limitations disallow programmatic instantiation of parameters
    p(ρ) = 5^2 * ρ
    β = 0.001
    tspan = (0.0, 1.0)
    dt = 0.25
    param_dt = 0.5
    params_tsteps = tspan[1]:param_dt:tspan[end]
    sys = create_system(g, t, dx, params)
    sys = structural_simplify(sys)
    prob = ODEProblem(
        sys,
        ones(length(sys.eqs)),
        tspan,
        Symbolics.scalarize(params .=> zeros(length(params))),
    )
    sol = solve(prob, Euler(), adaptive = false, dt = dt)

    loss(p) = sum(abs2, 2.0 .- solve(prob, Euler(), adaptive=false, dt=dt, p=p)[end]);
    η = 1e-1;

    #optimization loop
    p0 = ones(length(params));
    for i in 1:1000
        println(loss(p0))
        gs = gradient(loss, p0)[1]
        p0 -= η * gs;
    end
end
