@testset "lrsnash.jl" begin

    @testset "test 3 by 2 non-degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                           Player([1 0 0; 0 0 0]))
        NEs = [([1, 0, 0], [1, 0]),
               ([0, 1, 0], [0, 1]),
               ([0, 1, 0], [1//2, 1//2]),
               ([0, 0, 1], [1//2, 1//2])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

end
