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

end
