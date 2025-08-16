@testset "Testing Support Enumeration" begin

    function NEs_approx_equal(NEs1::Vector{NTuple{2,Vector{T1}}},
                              NEs2::Vector{NTuple{2,Vector{T2}}}) where {T1,T2}
        @test length(NEs1) == length(NEs2)
        @test T1 == T2
        for (actions1, actions2) in zip(NEs1, NEs2)
            for (action1, action2) in zip(actions1, actions2)
                @test action1 ≈ action2
            end
        end
    end

    @testset "test 3 by 2 non-degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([3.0 3.0; 2.0 5.0; 0.0 6.0]),
                           Player([3.0 2.0 3.0; 2.0 6.0 1.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([1.0 -1.0; -1.0 1.0; 0.0 0.0]),
                           Player([1.0 0.0 0.0; 0.0 0.0 0.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                           Player([1 0 0; 0 0 0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([1//1 -1//1; -1//1 1//1; 0//1 0//1]),
                           Player([1//1 0//1 0//1; 0//1 0//1 0//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([0//1, 1//1, 0//1], [0//1, 1//1])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

end
