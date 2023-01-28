@testset "homotopy_continuation.jl" begin

    @testset "3x2 game" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]

        NEs_computed = @inferred hc_solve(g, show_progress=false)
        @test isapprox_vecs_act_profs(NEs_computed, NEs)
    end

    @testset "2x2x2 game from McKelvey and McLennan" begin
        g = NormalFormGame((2, 2, 2))
        g[1, 1, 1] = 9, 8, 12
        g[2, 2, 1] = 9, 8, 2
        g[1, 2, 2] = 3, 4, 6
        g[2, 1, 2] = 3, 4, 4
        NEs = [
            ([1, 0], [1, 0], [1, 0]),
            ([0, 1], [0, 1], [1, 0]),
            ([1, 0], [0, 1], [0, 1]),
            ([0, 1], [1, 0], [0, 1]),
            ([0//1, 1//1], [1//3, 2//3], [1//3, 2//3]),
            ([1//4, 3//4], [1//1, 0//1], [1//4, 3//4]),
            ([1//2, 1//2], [1//2, 1//2], [1//1, 0//1]),
            ([1//4, 3//4], [1//2, 1//2], [1//3, 2//3]),
            ([1//2, 1//2], [1//3, 2//3], [1//4, 3//4])
        ]

        NEs_computed = @inferred hc_solve(g, show_progress=false)
        @test isapprox_vecs_act_profs(NEs_computed, NEs)

        ntofind = 1
        NEs_computed =
            @inferred hc_solve(g, ntofind=ntofind, show_progress=false)
        @test length(NEs_computed) == ntofind
        for i in 1:ntofind
            @test is_nash(g, NEs_computed[i])
        end
    end

    @testset "2x2x2 game from Nau, Canovas, and Hansen" begin
        payoff_profiles = [[3, 0, 2],
                           [0, 1, 0],
                           [0, 2, 0],
                           [1, 0, 0],
                           [1, 0, 0],
                           [0, 3, 0],
                           [0, 1, 0],
                           [2, 0, 3]]
        g = NormalFormGame(reshape(payoff_profiles, (2, 2, 2)))
        q = (-13 + sqrt(601)) / 24
        p = (9q - 1) / (7q + 2)
        r = (-3q + 2) / (q + 1)
        NEs = [([p, 1-p], [q, 1-q], [r, 1-r])]

        NEs_computed = @inferred hc_solve(g, show_progress=false)
        @test isapprox_vecs_act_profs(NEs_computed, NEs)
    end

    @testset "1-player game" begin
        g = NormalFormGame([[1], [2], [3]])
        @test_throws ArgumentError hc_solve(g)
        @test_throws ArgumentError hc_solve(g, ntofind=1)
    end

end
