@testset "bimatrix_generators.jl" begin

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

    @testset "unit_vector_game" begin
        k = 5
        g1 = @inferred(unit_vector_game(MersenneTwister(0), k; random=false))
        g2 = @inferred(unit_vector_game(MersenneTwister(0), k; random=false))

        for i in 1:2
            @test g1.players[i].payoff_array == g2.players[i].payoff_array
        end

        for i in 1:k
            n_0s = 0
            n_1s = 0
            for j in g1.players[1].payoff_array[:, i]
                if j == 0.0
                    n_0s += 1
                elseif j == 1.0
                    n_1s += 1
                end
            end
            @test n_0s == k - 1
            @test n_1s == 1
        end

        for k in 2:10
            @test pure_nash(unit_vector_game(k, random=false)) == []
        end

    end

end
