using Compat.LinearAlgebra
using Compat.Random
using Compat.Unicode
using Compat.DelimitedFiles
import Combinatorics: binomial

@testset "bimatrix_generators.jl" begin

    @testset "blotto_game" begin
        h = 3
        t = 7
        rho = 0.5

        g = @inferred(blotto_game(h, t, rho))

        @testset "test_nb_actions" begin
            nb_actions = (factorial(t+h-1)/
                          (factorial(h-1)*factorial(t)))

            @test g.nums_actions == (nb_actions, nb_actions)
        end

        @testset "test_constant_diagonal" begin
            for i=1:2
                d = diag(g.players[i].payoff_array)
                @test length(unique(d)) == 1
            end
        end

        @testset "test_seed" begin
            g1 = blotto_game(MersenneTwister(32), h, t, rho)
            g2 = blotto_game(MersenneTwister(32), h, t, rho)
            for i in 1:2
                @test g1.players[i].payoff_array == g2.players[i].payoff_array
            end
        end

    end

    @testset "ranking_game" begin
        g = @inferred(ranking_game(100))
        p1_array = g.players[1].payoff_array
        p2_array = g.players[2].payoff_array

        @testset "test_weakly_decreasing_rowise_payoffs" begin
            @test all(p1_array[:, 1:(end-1)] - p1_array[:, 2:end] .>=0)
            @test all(p2_array[:, 1:(end-1)] - p2_array[:, 2:end] .>=0)
        end

        @testset "test_elements_first_row" begin
            @test sum(g[1, 1]) == 1
            possible_elements = [0, 1, 0.5]
            @test all([value in possible_elements for value in p1_array[1, :]])
            @test all([value in possible_elements for value in p2_array[1, :]])
        end

        @testset "test_seed" begin
            g1 = ranking_game(MersenneTwister(32), 20)
            g2 = ranking_game(MersenneTwister(32), 20)
            for i in 1:2
                @test g1.players[i].payoff_array == g2.players[i].payoff_array
            end
        end

    end

    @testset "sgc_game" begin
        k = 2
        n = 4*k - 1
        g = @inferred(sgc_game(k))

        s = "
            0.750 0.750 1.000 0.500 0.500 1.000 0.000 0.500 0.000 0.500 0.000
            0.500 0.000 0.500 0.500 1.000 0.750 0.750 1.000 0.500 0.000 0.500
            0.000 0.500 0.000 0.500 0.000 0.500 1.000 0.500 0.500 1.000 0.750
            0.750 0.000 0.500 0.000 0.500 0.000 0.500 0.000 0.500 0.500 0.000
            0.500 0.000 0.500 0.000 0.750 0.000 0.000 0.750 0.000 0.000 0.000
            0.000 0.500 0.000 0.500 0.000 0.500 0.000 0.000 0.750 0.750 0.000
            0.000 0.000 0.000 0.000 0.500 0.000 0.500 0.000 0.500 0.000 0.000
            0.000 0.000 0.000 0.750 0.000 0.000 0.750 0.500 0.000 0.500 0.000
            0.500 0.000 0.000 0.000 0.000 0.000 0.000 0.750 0.750 0.000"
        payoffs = readdlm(IOBuffer(Unicode.normalize(s, stripcc=true)))
        payoffs = reshape(payoffs, (2, n^2))
        payoff_matrices = [reshape(payoffs[i, :], (n, n)) for i in 1:2]
        payoff_matrices[2] = transpose(payoff_matrices[2])

        for i in 1:2
            @test g.players[i].payoff_array == payoff_matrices[i]
        end
    end

    @testset "tournament_game" begin
        n, k = 6, 4
        m = binomial(n, k)
        g = @inferred(tournament_game(n, k))

        @testset "test_size" begin
            @test g.nums_actions == (n, m)
        end

        @testset "test_payoff_values" begin
            possible_values = [0., 1.]
            for player in g.players
                @test issubset(player.payoff_array, possible_values)
            end

            max_num_dominated_subsets = sum([binomial(i, k) for i in k:(n-1)])
            @test sum(g.players[1].payoff_array) <= max_num_dominated_subsets
            @test all(Compat.sum(g.players[2].payoff_array, dims=2) .== k)
        end

        @testset "test_seed" begin
            seed = 0
            g1 = tournament_game(n, k; seed=seed)
            g2 = tournament_game(n, k; seed=seed)

            for i in 1:2
                @test g1.players[i].payoff_array == g2.players[i].payoff_array
            end
        end

        @testset "test_throws_argument_error_too_large_inputs" begin
            n, k = 100, 50
            @test_throws ArgumentError tournament_game(n, k)
        end
    end

    @testset "unit_vector_game" begin
        n = 100
        g = @inferred(unit_vector_game(n))

        @testset "test_size" begin
            @test g.nums_actions == (n, n)
        end

        @testset "test_payoff_values" begin
            @test all(Compat.sum(g.players[1].payoff_array, dims=1) .== 1.)
        end

        @testset "test_avoid_pure_nash" begin
            NEs = pure_nash(unit_vector_game(n, avoid_pure_nash=true), tol=0.)
            @test length(NEs) == 0
        end

        @testset "test_seed" begin
            seed = 0
            n = 100
            g1 = unit_vector_game(MersenneTwister(seed), n)
            g2 = unit_vector_game(MersenneTwister(seed), n)
            for i in 1:2
                @test g1.players[i].payoff_array == g2.players[i].payoff_array
            end
        end

        @testset "test_redraw" begin
            seed = 6
            rng = MersenneTwister(seed)
            n = 2
            g = unit_vector_game(rng, n, avoid_pure_nash=true)
            NEs = pure_nash(g, tol=0.)
            @test length(NEs) == 0
        end

        @testset "test_throws_argument_error_avoid_pure_nash_n_1" begin
            n = 1
            @test_throws ArgumentError unit_vector_game(n, avoid_pure_nash=true)
        end
    end

end
