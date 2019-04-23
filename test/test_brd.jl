# ------------------------------- #
# Testing best response dynamics  #
# ------------------------------- #

using Random

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
        @test all(time_series(brd, ts_length) .<= N)
    end

    @testset "Testing KMR model" begin
        
        epsilon = 0.1
        seed = 1234
        kmr = KMR(payoff_matrix, N, epsilon)
        @test all(play(MersenneTwister(seed), kmr, init_action_dist,
                             num_reps=ts_length) .<= N)
        @test all(time_series(MersenneTwister(seed), kmr, ts_length,
                                    init_action_dist) .<= N)
        @test all(time_series(MersenneTwister(seed), kmr, ts_length) .<= N)
    end

    @testset "Testing sampling best response dynamics model" begin
        
        k = 2
        seed = 1234
        sbrd = SamplingBRD(payoff_matrix, N, k)
        @test all(play(MersenneTwister(seed), sbrd, init_action_dist,
                             num_reps=ts_length) .<= N)
        @test all(time_series(MersenneTwister(seed), sbrd, ts_length,
                                    init_action_dist) .<= N)
        @test all(time_series(MersenneTwister(seed), sbrd, ts_length) .<= N)
    end

    @testset "Testing argument errors" begin
        
        @testset "Non square payoff matrix" begin
            non_square_pay = ones((2, 3))

            @test_throws ArgumentError brd = BRD(non_square_pay, N)
            @test_throws ArgumentError kmr = KMR(non_square_pay, N, 0.1)
            @test_throws ArgumentError sbrd = SamplingBRD(non_square_pay, N, 2)
        end

        @testset "Invalid initial action distribution" begin
            invalid_action_dist_1 = [4, 0, 0]
            invalid_action_dist_2 = [3, 0]
            brd = BRD(payoff_matrix, N)

            @test_throws ArgumentError x = play(brd, invalid_action_dist_1)
            @test_throws ArgumentError x = play(brd, invalid_action_dist_2)
            @test_throws ArgumentError x = time_series(brd, 3, invalid_action_dist_1)
            @test_throws ArgumentError x = time_series(brd, 3, invalid_action_dist_2)
        end
    end
end