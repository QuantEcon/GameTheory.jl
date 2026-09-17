using Random

@testset "Testing Support Enumeration" begin

    function NEs_approx_equal(NEs1::Vector{NTuple{2,Vector{T1}}},
                              NEs2::Vector{NTuple{2,Vector{T2}}};
                              atol=0, rtol=1e-15) where {T1,T2}
        @test length(NEs1) == length(NEs2)
        @test T1 == T2
        for (actions1, actions2) in zip(NEs1, NEs2)
            for (action1, action2) in zip(actions1, actions2)
                @test isapprox(action1, action2, atol=atol, rtol=rtol)
            end
        end
    end

    @testset "test 3 by 2 non-degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([3.0 3.0; 2.0 5.0; 0.0 6.0]),
                           Player([3.0 2.0 3.0; 2.0 6.0 1.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(BigFloat)" begin
        g = NormalFormGame(Player([3.0 3.0; 2.0 5.0; 0.0 6.0]),
                           Player([3.0 2.0 3.0; 2.0 6.0 1.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.8, 0.2, 0.0], [2/3, 1/3]),
               ([0.0, 1/3, 2/3], [1/3, 2/3])]
        T = BigFloat
        g_BigFloat = NormalFormGame(T, g)
        NEs_BigFloat = [(T.(x), T.(y)) for (x, y) in NEs]
        NEs_computed = @inferred(support_enumeration(g_BigFloat))

        NEs_approx_equal(NEs_computed, NEs_BigFloat)
    end

    @testset "test 3 by 2 degenerate normal form game(Float)" begin
        g = NormalFormGame(Player([1.0 -1.0; -1.0 1.0; 0.0 0.0]),
                           Player([1.0 0.0 0.0; 0.0 0.0 0.0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                           Player([1 0 0; 0 0 0]))
        NEs = [([1.0, 0.0, 0.0], [1.0, 0.0]),
               ([0.0, 1.0, 0.0], [0.0, 1.0])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "test 3 by 2 degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([1//1 -1//1; -1//1 1//1; 0//1 0//1]),
                           Player([1//1 0//1 0//1; 0//1 0//1 0//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([0//1, 1//1, 0//1], [0//1, 1//1])]
        NEs_computed = @inferred(support_enumeration(g))

        NEs_approx_equal(NEs_computed, NEs)
    end

    @testset "N-player support enumeration" begin

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

            NEs_computed = @inferred support_enumeration(g)
            @test isapprox_vecs_act_profs(NEs_computed, NEs)

            # Explicit solver with options
            NEs_computed =
                @inferred support_enumeration(g, HCSolver(seed=UInt32(1234)))
            @test isapprox_vecs_act_profs(NEs_computed, NEs)

            ntofind = 1
            NEs_computed = @inferred support_enumeration(g, ntofind=ntofind)
            @test length(NEs_computed) == ntofind
            for i in 1:ntofind
                @test is_nash(g, NEs_computed[i])
            end

            # The per-support solver must be inferrable for the driver loop
            # to be; `@inferred support_enumeration` alone does not check
            # this, as the return type is fixed by the declaration of `NEs`
            supps = ntuple(i -> [1, 2], 3)
            sols = @inferred GameTheory._support_solutions(HCSolver(), g,
                                                          supps, [1, 2, 3])
            @test sols isa Vector{Vector{Float64}}
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

            NEs_computed = @inferred support_enumeration(g)
            @test isapprox_vecs_act_profs(NEs_computed, NEs)
        end

        @testset "3x2 games with explicit solver" begin
            # Non-degenerate game
            g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                               Player([3 2 3; 2 6 1]))
            NEs = support_enumeration(g)
            NEs_computed = @inferred support_enumeration(g, HCSolver())
            @test isapprox_vecs_act_profs(NEs_computed, NEs)

            # Degenerate game: player 2 is indifferent against action 1 of
            # player 1, so that ([1, 0, 0], [q, 1-q]) with 2/3 <= q <= 1 are
            # all Nash equilibria; only the pure one and the isolated mixed
            # one are found
            g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                               Player([3 2 3; 3 6 1]))
            NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
                   ([0//1, 1//3, 2//3], [1//3, 2//3])]
            NEs_computed = @inferred support_enumeration(g, HCSolver())
            @test isapprox_vecs_act_profs(NEs_computed, NEs)
        end

        @testset "degenerate game with all payoffs zero" begin
            g = NormalFormGame((2, 2, 2))
            NEs_computed = @inferred support_enumeration(g)
            @test length(NEs_computed) == 8
            for NE in NEs_computed
                @test all(x -> all(in((0, 1)), x), NE)
            end
        end

        @testset "random games" begin
            seed = 1234
            # Cross-check with hc_solve
            g = random_game(MersenneTwister(seed), (2, 2, 2, 2))
            NEs = hc_solve(g, show_progress=false, compile=false)
            NEs_computed = @inferred support_enumeration(g)
            @test isapprox_vecs_act_profs(NEs_computed, NEs)

            # Numbers of Nash equilibria verified against hc_solve
            for (nums_actions, num_NEs) in [((2, 2, 2, 2, 2), 3),
                                            ((3, 3, 2, 2), 5),
                                            ((3, 3, 3), 3)]
                g = random_game(MersenneTwister(seed), nums_actions)
                NEs_computed = @inferred support_enumeration(g)
                @test length(NEs_computed) == num_NEs
                for NE in NEs_computed
                    @test is_nash(g, NE)
                end
            end
        end

        @testset "1-player game" begin
            g = NormalFormGame([[1], [2], [3]])
            @test_throws ArgumentError support_enumeration(g)
            @test_throws ArgumentError support_enumeration(g, HCSolver())
        end

    end

end
