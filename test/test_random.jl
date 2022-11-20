using Random

@testset "Testing Random Games Generating" begin

    @testset "test random game" begin
        nums_actions = (2, 3, 4)
        g = @inferred random_game(nums_actions)
        @test g.nums_actions == nums_actions

        nums_actions = (4, 3)
        N = length(nums_actions)
        seed = 1234
        rngs = [MersenneTwister(seed) for i in 1:2]
        gs_Float = [random_game(rng, nums_actions) for rng in rngs]
        gs_Range = [random_game(rng, 0:10, nums_actions) for rng in rngs]
        for gs in [gs_Float, gs_Range]
            for i in 1:N
                @test gs[1].players[i].payoff_array ==
                      gs[2].players[i].payoff_array
            end
        end
    end

    @testset "test covariance game" begin
        nums_actions = (2, 3, 4)
        N = length(nums_actions)

        rho = 0.5
        g = covariance_game(nums_actions, rho)
        @test g.nums_actions == nums_actions

        rho = 1
        g = covariance_game(nums_actions, rho)
        for a in CartesianIndices(nums_actions)
            payoff_profile = g[a]
            for i in 1:(N-1)
                @test payoff_profile[i] ≈ payoff_profile[end]
            end
        end

        rho = -1/(N-1)
        g = covariance_game(nums_actions, rho)
        for a in CartesianIndices(nums_actions)
            payoff_profile = g[a]
            @test sum(payoff_profile) ≈ 0 atol=1e-10
        end

    end

    @testset "test random game value error" begin
        nums_actions = ()

        @test_throws ArgumentError random_game(nums_actions)

    end

    @testset "test covariance game value error" begin
        nums_actions = () #empty

        @test_throws ArgumentError covariance_game(nums_actions, 0.5)

        nums_actions = (2,) #length one

        @test_throws ArgumentError covariance_game(nums_actions, 0.5)

        nums_actions = (2, 3, 4)
        rho = 1.1 # > 1

        @test_throws ArgumentError covariance_game(nums_actions, rho)

        rho = -1. # < -1/(N-1)

        @test_throws ArgumentError covariance_game(nums_actions, rho)

    end

    @testset "random_pure_actions" begin
        nums_actions = (2, 3, 4)
        seed = 1234
        action_profiles =
            [random_pure_actions(MersenneTwister(seed), nums_actions)
             for i in 1:2]
        @test action_profiles[1] <= nums_actions
        @test action_profiles[2] == action_profiles[1]
    end

    @testset "random_mixed_actions" begin
        nums_actions = (2, 3, 4)
        seed = 1234
        action_profiles =
            [random_mixed_actions(MersenneTwister(seed), nums_actions)
             for i in 1:2]
        @test length.(action_profiles[1]) == nums_actions
        for i in 1:length(nums_actions)
            @test action_profiles[2][i] == action_profiles[1][i]
        end
    end

end
