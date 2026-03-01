using GameTheory: GAMPayoffVector

using Random

@testset "game_converters.jl" begin

    @testset "GAMPayoffVector" begin
        @testset "Golden: N=3" begin
            nums_actions = (2, 3, 4)
            N = length(nums_actions)
            na = prod(nums_actions)

            A1 = reshape(collect(1:na), nums_actions)
            A2 = reshape(collect(101:100+na), nums_actions)
            A3 = reshape(collect(201:200+na), nums_actions)

            payoffs1d = vcat(vec(A1), vec(A2), vec(A3))

            p = @inferred GAMPayoffVector(nums_actions, payoffs1d)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d
            @test num_players(p) == N

            payoffs4d = Array{Int,N+1}(undef, nums_actions..., N)
            payoffs4d[:, :, :, 1] .= A1
            payoffs4d[:, :, :, 2] .= A2
            payoffs4d[:, :, :, 3] .= A3

            g = NormalFormGame(payoffs4d)
            p = @inferred GAMPayoffVector(g)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d

            g_from_p = @inferred NormalFormGame(p)

            @test g_from_p.nums_actions == g.nums_actions
            for i in 1:N
                @test g_from_p.players[i].payoff_array ==
                      g.players[i].payoff_array
            end

            # Make an AbstractVector (SubArray) that equals payoffs1d
            payoffs1d_view = @view vcat([-999], payoffs1d, [999])[2:end-1]

            p = @inferred GAMPayoffVector(nums_actions, payoffs1d_view)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d
        end

        @testset "Round trip: N=$(length(ns))" for ns in [(4, 3), (2, 2, 3, 2)]
            N = length(ns)
            rng = MersenneTwister(12345)
            g = random_game(rng, 0:99, ns)
            p = @inferred GAMPayoffVector(g)
            g2 = @inferred NormalFormGame(p)

            p_BI = @inferred GAMPayoffVector(BigInt, g)
            g3 = @inferred NormalFormGame(Int, p_BI)

            for g_new in [g2, g3]
                @test g_new.nums_actions == g.nums_actions
                for i in 1:N
                    @test g_new.players[i].payoff_array ==
                          g.players[i].payoff_array
                end
            end
        end

        @testset "N=1" begin
            payoffs = [1., 2., 3.]
            nums_actions = (3,)

            p1 = GAMPayoffVector(nums_actions, payoffs)

            g = NormalFormGame(Player(payoffs))
            p2 = GAMPayoffVector(g)

            for p in [p1, p2]
                @test p.nums_actions == nums_actions
                @test p.payoffs == payoffs
            end
        end

        @testset "Invalid inputs" begin
            @test_throws ArgumentError GAMPayoffVector((2, 2), [1, 2, 3])
            @test_throws ArgumentError GAMPayoffVector((2, 0), Int[])
        end
    end

end
