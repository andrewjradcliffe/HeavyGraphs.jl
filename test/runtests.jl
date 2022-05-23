using HeavyGraphs
using Test

@testset "HeavyGraphs.jl" begin
    @testset "SimpleDiGraph: Constructors at 1 level of depth" begin
        x = SimpleDiGraph{UInt}()
        ks = [1:10;]
        get!.(() -> SimpleDiGraph{UInt}(), Ref(x), ks)
        @test keys(x) ⊆ ks
        @test all(values(x) .== Ref(SimpleDiGraph{UInt}()))
        @test length(x) == 10
        @test rlength(x) == 11
        @test depth(x) == 2
        @test size(x) == (1, 10)
        xs = map(x -> SimpleDiGraph{UInt}(), 1:5)
        setindex!.(Ref(x), xs, 11:15)
        @test length(x) == 15
        @test rlength(x) == 16
        @test depth(x) == 2
        @test size(x) == (1, 15)
        @test getindex.(Ref(x), 11:15) == xs
        x[10] = SimpleDiGraph()
        @test length(x) == 15
        @test rlength(x) == 16
        @test depth(x) == 2
        @test size(x) == (1, 15)
        @test get(() -> nothing, x, 5) isa SimpleDiGraph
        @test get(() -> nothing, x, 20) === nothing
        @test get!(() -> SimpleDiGraph{UInt}(), x, 16) isa SimpleDiGraph
        @test length(x) == 16
        @test rlength(x) == 17
        @test depth(x) == 2
        @test size(x) == (1, 16)
        @test haspath(x, 1)
        y = deepcopy(x)
        @test x == y
        @test setindex!(x, SimpleDiGraph(), 10) == y
    end
    @testset "SimpleDiGraph: Constructors at ≥ 2 levels of depth" begin
        x = SimpleDiGraph{UInt}()
        ks = [1:10;]
        get!(() -> SimpleDiGraph{UInt}(), x, ks...)
        @test length(x) == 1
        @test rlength(x) == 11
        @test depth(x) == 11
        @test size(x) == (1,1,1,1,1,1,1,1,1,1,1)
        @test size(x, 10) == 1
        @test haspath(x, 1,2,3,4,5,6,7,8,9,10)
        @test !haspath(x, 3,2,3,4,5)
        @test !haspath(x, 2,2,3,4,5,6,7,8,9,10)
        ks′ = [1:10;];
        ks′[1] = 2
        get!(() -> SimpleDiGraph{UInt}(), x, ks′...)
        @test length(x) == 2
        @test rlength(x) == 21
        @test depth(x) == 11
        @test size(x) == (1,2,2,2,2,2,2,2,2,2,2)
        @test size(x, 10) == 2
        @test haspath(x, 1,2,3,4,5,6,7,8,9,10)
        @test !haspath(x, 3,2,3,4,5)
        @test haspath(x, 2,2,3,4,5,6,7,8,9,10)
        y = deepcopy(x)
        @test x == y
        delete!(x, 1)
        @test length(x) == 1
        @test rlength(x) == 11
        @test depth(x) == 11
        @test size(x) == (1,1,1,1,1,1,1,1,1,1,1)
        @test !haspath(x, 1,2,3,4,5,6,7,8,9,10)
        @test haspath(x, 2,2,3,4,5,6,7,8,9,10)
    end
    @testset "SimpleDiGraph: pop!, push!, delete!, etc." begin
        x = SimpleDiGraph{UInt}()
        ks = [1:10;]
        get!(() -> SimpleDiGraph{UInt}(), x, ks...)
        ks′ = [1:10;];
        ks′[1] = 2
        get!(() -> SimpleDiGraph{UInt}(), x, ks′...)
        y = deepcopy(x)
        delete!(x, 1)
        @test length(x) == 1
        @test rlength(x) == 11
        @test depth(x) == 11
        @test size(x) == (1,1,1,1,1,1,1,1,1,1,1)
        @test size(x, 10) == 1
        @test !haspath(x, 1,2,3,4,5,6,7,8,9,10)
        @test haspath(x, 2,2,3,4,5,6,7,8,9,10)
        x′ = pop!(x)
        @test length(x) == 0
        @test rlength(x) == 1
        @test depth(x) == 1
        @test size(x) == (1,)
        @test isempty(x)
        push!(x, x′)
        @test length(x) == 1
        @test rlength(x) == 11
        @test depth(x) == 11
        @test size(x) == (1,1,1,1,1,1,1,1,1,1,1)
        @test !haspath(x, 1,2,3,4,5,6,7,8,9,10)
        @test haspath(x, 2,2,3,4,5,6,7,8,9,10)
    end
    @testset "SimpleDiGraph: growth" begin
        for T ∈ (Any, UInt, Int, Float64, Float32)
            # breadth
            d, b = 4, 50
            mat = reshape([1:200;], (d, b));
            es = Edges(IndexedEdge(i) for i = 1:d)
            x = grow(() -> SimpleDiGraph{T}(), es, eachcol(mat));
            @test length(x) == b
            @test rlength(x) == d * b + 1
            @test depth(x) == d + 1
            @test size(x) == (1,b,b,b,b)
            @test size(x, d + 1) == b
            @test maxbreadth(x) == b
            @test haspath(x, 1,2,3,4)
            @test haspath(x, 5,6,7,8)
            @test !haspath(x, 2,2,3,4)
            # depth
            d, b = 50, 4
            mat = reshape([1:200;], (d, b));
            es = Edges(IndexedEdge(i) for i = 1:d)
            x = grow(() -> SimpleDiGraph{T}(), es, eachcol(mat))
            @test length(x) == b
            @test rlength(x) == d * b + 1
            @test depth(x) == d + 1
            @test size(x) == (1,ntuple(_ -> b, Val(d))...)
            @test size(x, d + 1) == b
            @test maxbreadth(x) == b
            @test haspath(x, 1:50...)
            @test !haspath(x, 5,6,7,8)
            @test !haspath(x, 2,2,3,4)
        end
    end
    @testset "SimpleGraph: Constructors" begin
        #
        usg() = SimpleGraph{UInt}()
        x = SimpleGraph{UInt}()
        bget!(usg, x, 1,2)
        size(x)
        y = SimpleGraph{UInt}()
        y1 = bget!(usg, y, 1)
        y2 = bget!(usg, y1, 2)
        @test y == x
        @test x[1] == y[1]
        @test x[1, 2] == y[1, 2]
        #
        x = SimpleGraph()
        bget!(usg, x, 1)
        y = SimpleGraph()
        get!(usg, y, 1)
        @test x == y
        bget!(usg, x, 1,2)
        get!(usg, y, 1,2)
        @test y == x
        @test x[1] == y[1]
        @test x[1, 2] == y[1, 2]
        get!(usg, x, 1,2,3,4,5)
        x1 = x[1]
        @test x1.badj[1] == x
        x2 = x1[2]
        @test x2.badj[2] == x1
        x3 = x[1,2,3]
        @test x3 == x2[3]
        @test x[1,2,3,4] == x3[4]
        tmp_p = x
        tmp_n = x
        for i ∈ (1,2,3,4,5)
            tmp_n = tmp_p[i]
            @test tmp_n.badj[i] == tmp_p
            @test all(==(1), size(tmp_n))
            tmp_p = tmp_n
        end
        #
        get!(usg, x, 2,3,4,5,6)
        @test size(x) == (1, 2,2,2,2,2)
        tmp_p = x
        tmp_n = x
        for i ∈ (2,3,4,5,6)
            tmp_n = tmp_p[i]
            @test tmp_n.badj[i] == tmp_p
            @test all(==(1), size(tmp_n))
            tmp_p = tmp_n
        end
        x[3] = SimpleGraph()
        x3 = x[3]
        @test x3.badj[3] == x
        get!(usg, x, 4)
        x4 = x[4]
        @test x4.badj[4] == x
        @test setindex!(x, SimpleGraph(), 5) == x
    end
end
