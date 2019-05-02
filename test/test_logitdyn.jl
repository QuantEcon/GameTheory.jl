# ------------------------------------- #
# Testing Logit-Dynamics Response model #
# ------------------------------------- #

using Random

@testset "Testing logitdyn.jl" begin

    payoff_matrix = [4 0; 3 2]
    beta = 4.0
    g = NormalFormGame(payoff_matrix)
    lgdy = LogitDynamics(g, beta)

    seed = 1234
    init_actions = (1, 1)
    ts_length = 3
    @test all(play(MersenneTwister(seed), lgdy, init_actions) .<= 2)
    @test all(play(lgdy, init_actions) .<= 2)
    @test all(time_series(MersenneTwister(seed), lgdy, ts_length, init_actions) .<= 2)
    @test all(time_series(lgdy, ts_length, init_actions) .<= 2)
end 
