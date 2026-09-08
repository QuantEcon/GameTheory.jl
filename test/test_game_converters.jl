using GameTheory:
    PayoffVector, PlayerMajor, ProfileMajor, GAMPayoffVector, NFGPayoffVector,
    _player_block

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

    @testset "NFGPayoffVector" begin
        @testset "Golden: N=3" begin
            nums_actions = (2, 3, 4)
            N = length(nums_actions)
            na = prod(nums_actions)

            A1 = reshape(collect(1:na), nums_actions)
            A2 = reshape(collect(101:100+na), nums_actions)
            A3 = reshape(collect(201:200+na), nums_actions)

            # Profile-major: row-major vectorization of the (na, N) matrix
            payoffs1d = vec(permutedims(hcat(vec(A1), vec(A2), vec(A3))))

            p = @inferred NFGPayoffVector(nums_actions, payoffs1d)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d
            @test num_players(p) == N

            payoffs4d = Array{Int,N+1}(undef, nums_actions..., N)
            payoffs4d[:, :, :, 1] .= A1
            payoffs4d[:, :, :, 2] .= A2
            payoffs4d[:, :, :, 3] .= A3

            g = NormalFormGame(payoffs4d)
            p = @inferred NFGPayoffVector(g)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d

            g_from_p = @inferred NormalFormGame(p)

            @test g_from_p.nums_actions == g.nums_actions
            for i in 1:N
                @test g_from_p.players[i].payoff_array ==
                      g.players[i].payoff_array
            end

            payoffs1d_view = @view vcat([-999], payoffs1d, [999])[2:end-1]

            p = @inferred NFGPayoffVector(nums_actions, payoffs1d_view)

            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs1d
        end

        @testset "Round trip: N=$(length(ns))" for ns in [(4, 3), (2, 2, 3, 2)]
            N = length(ns)
            rng = MersenneTwister(12345)
            g = random_game(rng, 0:99, ns)
            p = @inferred NFGPayoffVector(g)
            g2 = @inferred NormalFormGame(p)

            p_BI = @inferred NFGPayoffVector(BigInt, g)
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

            p1 = NFGPayoffVector(nums_actions, payoffs)

            g = NormalFormGame(Player(payoffs))
            p2 = NFGPayoffVector(g)

            for p in [p1, p2]
                @test p.nums_actions == nums_actions
                @test p.payoffs == payoffs
            end
        end

        @testset "Invalid inputs" begin
            @test_throws ArgumentError NFGPayoffVector((2, 2), [1, 2, 3])
            @test_throws ArgumentError NFGPayoffVector((2, 0), Int[])
        end
    end

    @testset "Golden: N=2, both layouts" begin
        # 3x2 game with payoff profiles, in column-major order over
        # (a_1, a_2):
        #   (1,1): (3,2)  (2,1): (0,6)  (3,1): (2,1)
        #   (1,2): (1,3)  (2,2): (4,0)  (3,2): (5,4)
        nums_actions = (3, 2)
        g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]))

        # Profile-major: (payoffs at profile 1)..., (payoffs at profile 2)...
        payoffs_nfg = [3, 2, 0, 6, 2, 1, 1, 3, 4, 0, 5, 4]
        # Player-major: (payoffs to player 1)..., (payoffs to player 2)...
        payoffs_gam = [3, 0, 2, 1, 4, 5, 2, 6, 1, 3, 0, 4]

        @test NFGPayoffVector(g).payoffs == payoffs_nfg
        @test GAMPayoffVector(g).payoffs == payoffs_gam

        for (PV, payoffs) in [(NFGPayoffVector, payoffs_nfg),
                              (GAMPayoffVector, payoffs_gam)]
            p = PV(nums_actions, payoffs)
            g_from_p = NormalFormGame(p)
            for i in 1:2
                @test g_from_p.players[i].payoff_array ==
                      g.players[i].payoff_array
            end
        end
    end

    @testset "Layouts" begin
        @test GAMPayoffVector === PayoffVector{PlayerMajor}
        @test NFGPayoffVector === PayoffVector{ProfileMajor}
        @test_throws MethodError PayoffVector((3, 2), collect(1:12))

        nums_actions = (3, 2)
        p_gam = GAMPayoffVector(nums_actions, collect(1:12))
        p_nfg = NFGPayoffVector(nums_actions, [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12])

        @test p_gam isa GAMPayoffVector
        @test !(p_gam isa NFGPayoffVector)

        @testset "Conversion between layouts" begin
            p = @inferred NFGPayoffVector(p_gam)
            @test p.nums_actions == p_nfg.nums_actions
            @test p.payoffs == p_nfg.payoffs

            p = @inferred GAMPayoffVector(p_nfg)
            @test p.payoffs == p_gam.payoffs

            p = @inferred GAMPayoffVector(Float64, p_gam)
            @test p isa GAMPayoffVector{2,Float64}
            @test p.payoffs == p_gam.payoffs
            @test p.payoffs !== p_gam.payoffs
        end

        @testset "_player_block is a view" begin
            for p in [p_gam, p_nfg]
                b = @inferred _player_block(p, 2)
                @test size(b) == nums_actions
                @test Base.mightalias(b, p.payoffs)
                @test b == [7 10; 8 11; 9 12]

                g = NormalFormGame(p)
                @test !Base.mightalias(g.players[2].payoff_array, p.payoffs)
            end
        end
    end

end
