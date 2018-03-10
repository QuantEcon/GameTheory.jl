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

        @testset "test_payoff_generation" begin
            # WARNING: These results were computed by hand and depend on the
            # structure of the function being tested. Altering the structure of
            # the function will cause this test to fail.
            p1 = [0.16305 -0.353007 -0.353007 -0.353007
                  0.679107 0.16305 -0.353007 -0.353007
                  0.679107 0.679107 0.16305 -0.353007
                  0.679107 0.679107 0.679107 0.16305]
            p2 = [0.381845 -0.293291 -0.293291 -0.293291
                  1.05698 0.381845  -0.293291 -0.293291
                  1.05698 1.05698 0.381845 -0.293291
                  1.05698 1.05698 1.05698 0.381845]

            h = 2
            t = 3
            rho = 0.5
            rng = MersenneTwister(0)
            g = blotto_game(rng, h, t, rho)

            @test p1 ≈ g.players[1].payoff_array atol=1e-5
            @test p2 ≈ g.players[2].payoff_array atol=1e-5
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

        @testset "test_payoff_generation" begin
            # WARNING: These results were computed by hand and depend on the
            # structure of the function being tested. Altering the structure of
            # the function will cause this test to fail.
            g = ranking_game(MersenneTwister(0), 5)

            p1 = [0 0 0 0 0
                  0.88 -0.12 -0.12 -0.12 -0.12
                  0.72 0.72 0.72 -0.28 -0.28
                  0.66 0.66 0.66 0.66 -0.34
                  0.64 0.64 0.64 0.64 .64]

            p2 = [1 0 0 0 0
                  0.96 0.96 -0.04 -0.04 -0.04
                  0.8 0.8 -0.2 -0.2 -0.2
                  0.68 0.68 0.68 -0.32 -0.32
                  0.48 0.48 0.48 0.48 -0.52 ]

            @test p1 ≈ g.players[1].payoff_array
            @test p2 ≈ g.players[2].payoff_array
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
        payoffs = readdlm(IOBuffer(normalize_string(s, stripcc=true)))
        payoffs = reshape(payoffs, (2, n^2))
        payoff_matrices = [reshape(payoffs[i, :], (n, n)) for i in 1:2]
        payoff_matrices[2] = transpose(payoff_matrices[2])

        for i in 1:2
            @test g.players[i].payoff_array == payoff_matrices[i]
        end
    end

    @testset "tournament_game" begin
        n, k = 6, 4
        seed = 0
        g1 = @inferred(tournament_game(n, k; seed=seed))
        g2 = @inferred(tournament_game(n, k; seed=seed))

        for i in 1:2
            @test g1.players[i].payoff_array == g2.players[i].payoff_array
        end

        @test g1.nums_actions == (n, binomial(n, k))

        @test all(sum(g1.players[2].payoff_array, 2) .== k)

    end

    @testset "unit_vector_game" begin
        n = 100
        g = @inferred(unit_vector_game(n))

        @testset "test_size" begin
            @test g.nums_actions == (n, n)
        end

        @testset "test_payoff_values" begin
            @test all(sum(g.players[1].payoff_array, 1) .== 1.)
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

        @testset "test_raises_value_error_avoid_pure_nash_n_1" begin
            n = 1
            @test_throws ArgumentError unit_vector_game(n, avoid_pure_nash=true)
        end
    end

end
