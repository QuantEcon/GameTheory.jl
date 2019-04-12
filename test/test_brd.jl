# ------------------------------- #
# Testing best response dynamics  #
# ------------------------------- #


@testset "Testing brd.jl" begin
    
    payoff_matrix = [4 0; 3 2]
    N = 4
    ts_length = 3
    init_action_dist = [4, 0]

    @testset "Testing best response dynamics model" begin
        
        brd = BRD(payoff_matrix, N)
        @test @inferred(play(brd, init_action_dist, num_reps=ts_length)) ==
              [4, 0]
        @test @inferred(time_series(brd, ts_length, init_action_dist)) ==
              [4 4 4; 0 0 0]
    end

    @testset "Testing KMR model" begin
        
        epsilon = 0.1
        kmr = KMR(payoff_matrix, N, epsilon)
        action_dist = play(kmr, init_action_dist, num_reps=ts_length)
        @test sum(action_dist) == 4
        series = time_series(kmr, ts_length, init_action_dist)
        for t in 1:3
            @test sum(series[:, t]) == 4
        end
    end

    @testset "Testing sampling best response dynamics model" begin
        
        k = 2
        sbrd = SamplingBRD(payoff_matrix, N, k)
        action_dist = play(sbrd, init_action_dist, num_reps=ts_length)
        @test sum(action_dist) == 4
        series = time_series(sbrd, ts_length, init_action_dist)
        for t in 1:3
            @test sum(series[:, t]) == 4
        end
    end
end