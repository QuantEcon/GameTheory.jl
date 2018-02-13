@testset "Testing Support Enumeration" begin

    @testset "test 3 by 2 non-degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([3.0 3.0; 2.0 5.0; 0.0 6.0]),
                           Player([3.0 2.0 3.0; 2.0 6.0 1.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: AbstractFloat
            end
        end
    end

    @testset "test 3 by 2 non-degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: AbstractFloat
            end
        end
    end

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: Rational
            end
        end
    end

    @testset "test 3 by 2 degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([1.0 -1.0; -1.0 1.0; 0.0 0.0]),
                           Player([1.0 0.0 0.0; 0.0 0.0 0.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: AbstractFloat
            end
        end
    end

    @testset "test 3 by 2 degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                           Player([1 0 0; 0 0 0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: AbstractFloat
            end
        end
    end

    @testset "test 3 by 2 degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([1//1 -1//1; -1//1 1//1; 0//1 0//1]),
                           Player([1//1 0//1 0//1; 0//1 0//1 0//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([0//1, 1//1, 0//1], [0//1, 1//1])]

        for (actions_computed, actions) in zip(support_enumeration(g), NEs)
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
                @test eltype(action_computed) <: Rational
            end
        end
    end


end
