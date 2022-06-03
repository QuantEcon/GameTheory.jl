@testset "Testing nondegenerate bimatrix games" begin
    @testset "Test 2 by 2 normal form game with 3 equilibria" begin
        g = NormalFormGame(Player([1 0; 0 1]), Player([1 0; 0 1]))
        NEs = [([1.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0], [0.0, 1.0]),
               ([0.5, 0.5], [0.5, 0.5])]
        NEs_computed = @inferred(vertex_enumeration(g))
        @test sort(NEs) == sort(NEs_computed)
    end

    @testset "Test 2 by 3 normal form game with 1 equilibrium" begin
        g = NormalFormGame(Player([5 4 3; 4 1 2]), Player([5 1; 4 3; 3 2]))
        NEs = [([1.0, 0.0], [1.0, 0.0, 0.0])]
        NEs_computed = @inferred(vertex_enumeration(g))
        @test NEs == NEs_computed
    end

    @testset "Test Matching Pennies" begin
        MP = [1 -1; -1 1]
        g = NormalFormGame(Player(MP), Player(-MP))
        NEs = [([0.5, 0.5], [0.5, 0.5])]
        NEs_computed = @inferred(vertex_enumeration(g))
        @test NEs == NEs_computed
    end
    @testset "Test 3 by 2 normal form game where polyhedron features imprecision error" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]), Player([3 2 3; 2 6 1]))
        NEs = sort([([1.0, 0.0, 0.0], [1.0, 0.0]),
                    ([0.8, 0.2, 0.0], [2/3, 1/3]),
                    ([0.0, 1/3, 2/3], [1/3, 2/3])])
        NEs_computed = sort(@inferred(vertex_enumeration(g)))
        NExs = []
        NEys = []
        for (x,y) in NEs
            push!(NExs, x)
            push!(NEys, y)
        end
        NExs_computed = []
        NEys_computed = []
        for (x,y) in NEs_computed
            push!(NExs_computed, x)
            push!(NEys_computed, y)
        end
        @test NExs ≈ NExs_computed && NEys ≈ NEys_computed
    end
end


import GameTheory: is_nondegenerate

@testset "Testing degenerate game errors" begin
    @testset "Test degenerate 2 by 2 normal form game throws error" begin
        A = [1 1; 1 1]
        g = NormalFormGame(Player(A), Player(A))
        @test is_nondegenerate(g) == false
        @test_throws ErrorException("The vertex enumeration algorithm will not yield a solution for degenerate games.") vertex_enumeration(g)
    end
    @testset "Test that 3 player game throws error" begin
        g = random_game((2, 2, 2))
        @test_throws AssertionError vertex_enumeration(g)
        @test_throws ErrorException("Not a bimatrix game.") is_nondegenerate(g)
    end
end
