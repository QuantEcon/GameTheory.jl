@testset "Testing Support Enumeration" begin

    @testset "test 3 by 2 non-degenerate normal form game" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                            Player([3 2 3; 2 6 1]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
                ([0.8, 0.2, 0.0], [2/3, 1/3]),
                ([0.0, 1/3, 2/3], [1/3, 2/3])]

        for (actions_computed, actions) in zip(NEs, support_enumeration(g))
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
            end
        end

    end

    @testset "test 3 by 2 degenerate normal form game" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                            Player([1 0 0; 0 0 0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
                ([0.0, 1.0, 0.0], [0.0, 1.0])]

        for (actions_computed, actions) in zip(NEs, support_enumeration(g))
            for (action_computed, action) in zip(actions_computed, actions)
                @test action_computed ≈ action
            end
        end
    end


end
