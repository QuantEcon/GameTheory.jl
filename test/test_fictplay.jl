# ------------------------- #
# Testing fictitious play   #
# ------------------------- #


using Distributions
using Random


@testset "Testing fictplay.jl" begin

    matching_pennies_bimatrix = Array{Float64}(undef, 2, 2, 2)
    matching_pennies_bimatrix[:, 1, 1] = [1, -1]
    matching_pennies_bimatrix[:, 1, 2] = [-1, 1]
    matching_pennies_bimatrix[:, 2, 1] = [-1, 1]
    matching_pennies_bimatrix[:, 2, 2] = [1, -1]
    g = NormalFormGame(matching_pennies_bimatrix)

    gain = 0.1
    init_actions = (1,1)
    init_actions_mixed = ([1.0, 0.0], [1.0, 0.0])

    function vector_approximate_equal(vec1::NTuple{2,Vector{T1}},
                                      vec2::NTuple{2,Vector{T2}}) where {T1,T2}
        @test T1 == T2
        for (x1, x2) in zip(vec1, vec2)
            @test length(x1) == length(x2)
            for (xx1, xx2) in zip(x1, x2)
                @test xx1 ≈ xx2
            end
        end
    end

    function matrix_approximate_equal(mat1::NTuple{2,Matrix{T1}},
                                      mat2::NTuple{2,Matrix{T2}}) where {T1,T2}
        @test T1 == T2
        for (x1, x2) in zip(mat1, mat2)
            @test size(x1) == size(x2)
            row, col = size(x1)
            for i in 1:row
                for j in 1:col
                    @test x1[i, j] ≈ x2[i, j]
                end
            end
        end
    end

    @testset "AbstractFictitiousPlay from AbstractFictitiousPlay" begin
        fp = FictitiousPlay(g)
        fp_g = FictitiousPlay(fp)

        @test fp_g.players == fp.players
        @test fp_g.nums_actions == fp.nums_actions
        @test fp_g.gain == fp.gain

        normal = Normal()
        sfp = StochasticFictitiousPlay(g, normal)
        sfp_g1 = StochasticFictitiousPlay(sfp, normal)
        sfp_g2 = StochasticFictitiousPlay(sfp)

        @test sfp_g1.players == sfp.players
        @test sfp_g2.players == sfp.players
        @test sfp_g1.nums_actions == sfp.nums_actions
        @test sfp_g2.nums_actions == sfp.nums_actions
        @test sfp_g2.d == sfp.d
        @test sfp_g1.gain == sfp.gain
        @test sfp_g2.gain == sfp.gain
    end

    @testset "Testing fictitious play" begin

        fp_dec = FictitiousPlay(g)
        x = play(fp_dec, init_actions)
        x_mixed = play(fp_dec, init_actions_mixed)
        x_series = time_series(fp_dec, 3, init_actions)
        x_series_mixed = time_series(fp_dec, 3, init_actions_mixed)
        x_des = ([1.0, 0.0], [0.5, 0.5])
        x_series_des = ([1.0 1.0 1.0; 0.0 0.0 0.0], [1.0 0.5 1/3; 0.0 0.5 2/3])

        fp_con = FictitiousPlay(g, ConstantGain(gain))
        y = play(fp_con, init_actions)
        y_mixed = play(fp_con, init_actions_mixed)
        y_series = time_series(fp_con, 3, init_actions)
        y_series_mixed = time_series(fp_con, 3, init_actions_mixed)
        y_des = ([1.0, 0.0], [0.9, 0.1])
        y_series_des = (
                        [1.0 1.0 1.0; 0.0 0.0 0.0], [1.0 0.9 0.81; 0.0 0.1 0.19]
                       )

        vector_approximate_equal(x, x_des)
        vector_approximate_equal(x_mixed, x_des)
        vector_approximate_equal(y, y_des)
        vector_approximate_equal(y_mixed, y_des)
        matrix_approximate_equal(x_series, x_series_des)
        matrix_approximate_equal(x_series_mixed, x_series_des)
        matrix_approximate_equal(y_series, y_series_des)
        matrix_approximate_equal(y_series_mixed, y_series_des)

        seed = 1234
        for i in 1:2
            @test all(play(MersenneTwister(seed), fp_dec)[i] .<= 1.0)
            @test all(play(fp_dec)[i] .<= 1.0)
            @test all(play(MersenneTwister(seed), fp_con)[i] .<= 1.0)
            @test all(play(fp_con)[i] .<= 1.0)
            @test all(time_series(MersenneTwister(seed), fp_dec, 3)[i] .<= 1.0)
            @test all(time_series(fp_dec, 3)[i] .<= 1.0)
            @test all(time_series(MersenneTwister(seed), fp_con, 3)[i] .<= 1.0)
            @test all(time_series(fp_con, 3)[i] .<= 1.0)
        end
    end

    @testset "Testing stochastic fictitious play" begin

        normal = Normal()  #standard normal distribution
        seed = 1234

        sfp_dec = StochasticFictitiousPlay(g, normal)
        x = [play(MersenneTwister(seed), sfp_dec, init_actions) for i in 1:2]
        x_mixed = [play(MersenneTwister(seed), sfp_dec, init_actions_mixed) for i in 1:2]
        x_series = [time_series(MersenneTwister(seed), sfp_dec, 3, init_actions)
                    for i in 1:2]
        x_series_mixed = [time_series(MersenneTwister(seed), sfp_dec, 3,
                          init_actions_mixed) for i in 1:2]

        sfp_con = StochasticFictitiousPlay(g, normal, ConstantGain(gain))
        y = [play(MersenneTwister(seed), sfp_con, init_actions) for i in 1:2]
        y_mixed = [play(MersenneTwister(seed), sfp_con, init_actions_mixed) for i in 1:2]
        y_series = [time_series(MersenneTwister(seed), sfp_con, 3, init_actions)
                    for i in 1:2]
        y_series_mixed = [time_series(MersenneTwister(seed), sfp_con, 3,
                          init_actions_mixed) for i in 1:2]

        vector_approximate_equal(x[1], x[2])
        vector_approximate_equal(x_mixed[1], x_mixed[2])
        vector_approximate_equal(y[1], y[2])
        vector_approximate_equal(y_mixed[1], y_mixed[2])
        matrix_approximate_equal(x_series[1], x_series[2])
        matrix_approximate_equal(x_series_mixed[1], x_series_mixed[2])
        matrix_approximate_equal(y_series[1], y_series[2])
        matrix_approximate_equal(y_series_mixed[1], y_series_mixed[2])

        for i in 1:2
            @test all(play(MersenneTwister(seed), sfp_dec)[i] .<= 1.0)
            @test all(play(MersenneTwister(seed), sfp_con)[i] .<= 1.0)
            @test all(time_series(MersenneTwister(seed), sfp_dec, 3)[i] .<= 1.0)
            @test all(time_series(MersenneTwister(seed), sfp_con, 3)[i] .<= 1.0)
        end
    end
end