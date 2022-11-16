import GameTheory: is_nondegenerate

@testset "Testing Vertex Enumeration" begin

    @testset "Test 2 by 2 normal form game with 3 equilibria" begin
        g = NormalFormGame(Player([1 0; 0 1]), Player([1 0; 0 1]))
        NEs = [([1.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0], [0.0, 1.0]),
               ([0.5, 0.5], [0.5, 0.5])]
        NEs_computed = @inferred(vertex_enumeration(g))
        NEs_approx_equal(sort(NEs), sort(NEs_computed))
    end


    @testset "Test 2 by 3 normal form game with 1 equilibrium" begin
        g = NormalFormGame(Player([5 4 3; 4 1 2]), Player([5 1; 4 3; 3 2]))
        NEs = [([1.0, 0.0], [1.0, 0.0, 0.0])]
        NEs_computed = @inferred(vertex_enumeration(g))
        NEs_approx_equal(sort(NEs), sort(NEs_computed))
    end


    @testset "Test Matching Pennies" begin
        MP = [1 -1; -1 1]
        g = NormalFormGame(Player(MP), Player(-MP))
        NEs = [([0.5, 0.5], [0.5, 0.5])]
        NEs_computed = @inferred(vertex_enumeration(g))
        NEs_approx_equal(NEs, NEs_computed)
    end


    @testset "Test Matching Pennies (rational)" begin
        MP = [1//1 -1//1; -1//1 1//1]
        g = NormalFormGame(Player(MP), Player(-MP))
        NEs = [([1//2, 1//2], [1//2, 1//2])]
        NEs_computed = @inferred(vertex_enumeration(g))
        NEs_approx_equal(NEs, NEs_computed)
    end
    

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(vertex_enumeration(g))
        NEs_approx_equal(NEs, NEs_computed)
    end

    @testset "Test 3 by 2 game where polyhedron has imprecision error" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]), Player([3 2 3; 2 6 1]))
        NEs = sort([([1.0, 0.0, 0.0], [1.0, 0.0]),
                    ([0.8, 0.2, 0.0], [2/3, 1/3]),
                    ([0.0, 1/3, 2/3], [1/3, 2/3])])
        NEs_computed = sort(@inferred(vertex_enumeration(g)))
        NEs_approx_equal(NEs, NEs_computed)
    end

    
    @testset "Test degenerate 2 by 2 normal form game throws error" begin
        A = [1 1; 1 1]
        g = NormalFormGame(Player(A), Player(A))
        @test is_nondegenerate(g) == false
        @test_throws ErrorException("The vertex enumeration algorithm will"*
        " not yield a solution for degenerate games.") vertex_enumeration(g)
    end

    
    @testset "Test that 3 player game throws error" begin
        g = random_game((2, 2, 2))
        @test_throws ErrorException("Not a bimatrix game.") is_nondegenerate(g)
    end
end
