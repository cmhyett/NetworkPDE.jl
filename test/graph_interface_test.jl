@testset "Component construction" begin
    # empty metadata
    e1 = Edge(1);
    @test e1 isa Edge
    @test e1.id isa Int
    @test e1.metadata isa Dict
    
    v1 = Vertex(1);
    @test v1 isa Vertex
    @test v1.id isa Int
    @test v1.metadata isa Dict

    # default constructor
    e1 = Edge(1, Dict(:L=>2))
    @test e1 isa Edge
    @test e1.id isa Int
    @test e1.metadata isa Dict

    v1 = Vertex(1, Dict(:is_flux => true));
    @test v1 isa Vertex
    @test v1.id isa Int
    @test v1.metadata isa Dict
end

@testset "Network construction" begin
    # empty construction
    @test Network() isa Network
    
    v1 = Vertex(1, Dict(:is_flux => true, :flux => zeros(3)));
    v2 = Vertex(2, Dict(:is_flux => true, :flux => zeros(3)));
    v3 = Vertex(3, Dict(:is_flux => true, :flux => zeros(3)));
    e1 = Edge(1, Dict(:L => 100.0, :λ => 0.11, :D => 0.25));
    e2 = Edge(2, Dict(:L => 50.0, :λ => 0.11, :D => 0.25));
    g = Network([e1, e2], [n1, n2, n3])
    @test g isa Network
end

@testset "Network Utilities" begin
    v1 = Vertex(1, Dict(:is_flux => true, :flux => zeros(3)));
    v2 = Vertex(2, Dict(:is_flux => true, :flux => zeros(3)));
    v3 = Vertex(3, Dict(:is_flux => true, :flux => zeros(3)));
    e1 = Edge(1, Dict(:L => 100.0, :λ => 0.11, :D => 0.25));
    e2 = Edge(2, Dict(:L => 50.0, :λ => 0.11, :D => 0.25));
    g = Network([e1, e2], [n1, n2, n3])
    @test g isa Network
end

