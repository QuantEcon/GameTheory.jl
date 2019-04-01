# ------------------------- #
# Testing fictitious play   #
# ------------------------- #


using Distributions


@testset "Testing fictplay.jl" begin

    matching_pennies_bimatrix = Array{Float64}(undef, 2, 2, 2)
    matching_pennies_bimatrix[:, 1, 1] = [1, -1]
    matching_pennies_bimatrix[:, 1, 2] = [-1, 1]
    matching_pennies_bimatrix[:, 2, 1] = [-1, 1]
    matching_pennies_bimatrix[:, 2, 2] = [1, -1]
    g = NormalFormGame(matching_pennies_bimatrix)

    gain = 0.1
    init_actions = (1,1)

    @testset "Testing fictitious play" begin

        fp_dec = FictitiousPlay(g)
        x = play(fp_dec, init_actions)
        x_series = time_series(fp_dec, 3, init_actions)
        x_des = ([1.0, 0.0], [0.5, 0.5])
        x_series_des = ([1.0 1.0 1.0; 0.0 0.0 0.0], [1.0 0.5 1/3; 0.0 0.5 2/3])

        fp_con = FictitiousPlay(g, ConstantGain(gain))
        y = play(fp_con, init_actions)
        y_series = time_series(fp_con, 3, init_actions)
        y_des = ([1.0, 0.0], [0.9, 0.1])
        y_series_des = (
                        [1.0 1.0 1.0; 0.0 0.0 0.0], [1.0 0.9 0.81; 0.0 0.1 0.19]
                       )

        for (i, j) in ((1, 1), (1, 2), (2, 1), (2, 2))
            @test x[i][j] ≈ x_des[i][j]
            @test y[i][j] ≈ y_des[i][j]
        end

        for k in 1:2
            for i in 1:2
                for j in 1:3
                    @test x_series[k][i, j] ≈ x_series_des[k][i, j]
                    @test y_series[k][i, j] ≈ y_series_des[k][i, j]
                end
            end
        end
    end

    @testset "Testing stochastic fictitious play" begin

        normal = Normal()  #standard normal distribution

        sfp_dec = StochasticFictitiousPlay(g, normal)
        sfp_con = StochasticFictitiousPlay(g, normal, ConstantGain(gain))

        x_dec = play(sfp_dec, init_actions)
        @test sum(x_dec[1]) ≈ 1
        @test sum(x_dec[2]) ≈ 1
        x_con = play(sfp_con, init_actions)
        @test sum(x_con[1]) ≈ 1
        @test sum(x_con[2]) ≈ 1

        y_dec = time_series(sfp_dec, 3, init_actions)
        for t in 1:3
            @test sum(y_dec[1][:,t]) ≈ 1
            @test sum(y_dec[2][:,t]) ≈ 1
        end
        y_con = time_series(sfp_con, 3, init_actions)
        for t in 1:3
            @test sum(y_con[1][:,t]) ≈ 1
            @test sum(y_con[2][:,t]) ≈ 1
        end
    end
end