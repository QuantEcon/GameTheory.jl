using GameTheory: BestResponsePolytope
using Polyhedra
using CDDLib
using LRSLib

@testset "Testing Vertex Enumeration" begin

    @testset "3x2 nondegenerate game" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        g_F = NormalFormGame(Float64, g)
        g_BF = NormalFormGame(BigFloat, g)
        g_Rat = NormalFormGame(Rational{BigInt}, g)
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]

        NEs_computed = @inferred vertex_enumeration(g)
        @test isapprox_vecs_act_profs(NEs_computed, NEs)

        for h in [g, g_F, g_BF]
            for plib in [default_library(2, Float64), CDDLib.Library()]
                NEs_computed = @inferred vertex_enumeration(h, plib=plib)
                @test isapprox_vecs_act_profs(NEs_computed, NEs)
            end
        end

        for plib in [default_library(2, Rational{BigInt}),
                     CDDLib.Library(:exact), LRSLib.Library()]
            NEs_computed = @inferred vertex_enumeration(g_Rat, plib=plib)
            @test isapprox_vecs_act_profs(NEs_computed, NEs)
        end
    end

    @testset "3x2 degenerate game" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 3 6 1]))

        # Not guaranteed to find all equilibria
        NEs_computed = @inferred vertex_enumeration(g)
        for NE in NEs_computed
            @test is_nash(g, NE)
        end
    end

    @testset "BestResponsePolytope" begin
        # From von Stengel 2007 in Algorithmic Game Theory
        A = [3 3; 2 5; 0 6]
        B = [3 2 3; 2 6 1]
        g = NormalFormGame(Player(A), Player(B))

        # Expected vertices keyed by (sorted) labelings
        expected = Dict{NTuple{3,Int}, Vector{Float64}}(
            (1, 2, 3) => [0.0, 0.0, 0.0],           # 0
            (2, 3, 4) => [1/3, 0.0, 0.0],           # a
            (3, 4, 5) => [2/7, 1/14, 0.0],          # b
            (1, 3, 5) => [0.0, 1/6, 0.0],           # c
            (1, 4, 5) => [0.0, 1/8, 1/4],           # d
            (1, 2, 4) => [0.0, 0.0, 1/3],           # e
        )

        i = 1
        brp = BestResponsePolytope(g.players[3-i];
                                   idx=i, plib=CDDLib.Library())
        vertices = points(vrep(brp.poly))

        @test length(vertices) == length(expected)

        for (idx, pt) in zip(eachindex(vertices), vertices)
            indices = incidenthalfspaceindices(brp.poly, idx)  # Vector of Polyhedra.Index
            labeling = sort(ntuple(i -> indices[i].value, 3))
            @test haskey(expected, labeling)
            @test isapprox(pt, expected[labeling])
        end
    end

end
