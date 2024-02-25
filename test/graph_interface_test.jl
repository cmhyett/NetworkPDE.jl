import NetworkPDE as npde

@testset "Component construction" begin
    # empty metadata
    e1 = npde.Edge(1)
    @test e1 isa npde.Edge
    @test e1.id isa Int
    @test e1.metadata isa Dict

    v1 = npde.Vertex(1)
    @test v1 isa npde.Vertex
    @test v1.id isa Int
    @test v1.metadata isa Dict

    # default constructor
    e1 = npde.Edge(1, Dict(:L => 2))
    @test e1 isa npde.Edge
    @test e1.id isa Int
    @test e1.metadata isa Dict

    v1 = npde.Vertex(1, Dict(:is_flux => true))
    @test v1 isa npde.Vertex
    @test v1.id isa Int
    @test v1.metadata isa Dict
end

@testset "Network construction" begin
    # empty construction
    @test npde.Network() isa npde.Network

    v1 = npde.Vertex(1, Dict(:is_flux => true, :flux => zeros(3)))
    v2 = npde.Vertex(2, Dict(:is_flux => true, :flux => zeros(3)))
    v3 = npde.Vertex(3, Dict(:is_flux => true, :flux => zeros(3)))
    e1 = npde.Edge(1, Dict(:L => 100.0, :λ => 0.11, :D => 0.25))
    e2 = npde.Edge(2, Dict(:L => 50.0, :λ => 0.11, :D => 0.25))
    g = npde.Network([e1, e2], [v1, v2, v3])
    @test g isa npde.Network
end

@testset "Network Utilities" begin
    v1 = npde.Vertex(1, Dict(:is_flux => true, :flux => zeros(3)))
    v2 = npde.Vertex(2, Dict(:is_flux => true, :flux => zeros(3)))
    v3 = npde.Vertex(3, Dict(:is_flux => true, :flux => zeros(3)))
    e1 = npde.Edge(1, Dict(:L => 100.0, :λ => 0.11, :D => 0.25))
    e2 = npde.Edge(2, Dict(:L => 50.0, :λ => 0.11, :D => 0.25))
    g = npde.Network([e1, e2], [v1, v2, v3])
    @test g isa npde.Network

    @test npde.get_metadata(v1, :is_flux) == true
    @test npde.get_metadata(e1, :L) == 100.0

    npde.add_metadata(v1, :A, 3.14)
    @test npde.get_metadata(v1, :A) == 3.14
    npde.add_metadata(v1, Dict(:B => 6.28))
    @test npde.get_metadata(v1, :B) == 6.28

    @test npde.get_edge(g, 1) isa npde.Edge
    @test npde.get_vertex(g, 2) isa npde.Vertex

    nedges = npde.num_edges(g)
    nvertices = npde.num_vertices(g)
    @test nedges == 2
    @test nvertices == 3
    npde.add_edge(g, npde.Edge(3))
    @test npde.num_edges(g) == 3
end
