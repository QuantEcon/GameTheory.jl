using GameTheory:
    PayoffVector, PlayerMajor, ProfileMajor, GAMPayoffVector, NFGPayoffVector,
    _player_block

using Random

# Element type that the writers do not accept
struct UnsupportedReal <: Real end

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

    @testset "read_gam/write_gam" begin
        same_game(g1, g2) =
            g1.nums_actions == g2.nums_actions &&
            all(g1.players[i].payoff_array == g2.players[i].payoff_array
                for i in 1:num_players(g1))

        @testset "Golden: N=2" begin
            g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]))
            s = "2\n3 2\n\n3 0 2 1 4 5 2 6 1 3 0 4\n"

            @test gam_string(g) == s
            @test sprint(write_gam, g) == s
            @test same_game(parse_gam(s), g)
            @test same_game(read_gam(IOBuffer(s)), g)

            p = GAMPayoffVector((3, 2), [3, 0, 2, 1, 4, 5, 2, 6, 1, 3, 0, 4])
            @test gam_string(p) == s
        end

        @testset "Golden: N=3" begin
            nums_actions = (2, 2, 2)
            g = NormalFormGame(Int, nums_actions)
            for (a, payoffs) in [((1, 1, 1), (0, 8, 16)), ((2, 1, 1), (1, 9, 17)),
                                 ((1, 2, 1), (2, 10, 18)), ((2, 2, 1), (3, 11, 19)),
                                 ((1, 1, 2), (4, 12, 20)), ((2, 1, 2), (5, 13, 21)),
                                 ((1, 2, 2), (6, 14, 22)), ((2, 2, 2), (7, 15, 23))]
                g[a...] = payoffs
            end
            s = "3\n2 2 2\n\n" * join(0:23, ' ') * '\n'

            @test gam_string(g) == s
            @test same_game(parse_gam(s), g)
        end

        @testset "Round trip: N=$(length(ns)), $S" for ns in [(4, 3), (2, 2, 3, 2)],
                                                        S in [0:99, Float64]
            rng = MersenneTwister(12345)
            g = random_game(rng, S, ns)

            @test same_game(parse_gam(gam_string(g)), g)

            mktempdir() do dir
                path = joinpath(dir, "game.gam")
                @test write_gam(path, g) === nothing
                @test read(path, String) == gam_string(g)
                @test same_game(read_gam(path), g)
                @test same_game(read_gam(eltype(g.players[1].payoff_array), path), g)
            end
        end

        @testset "Element type" begin
            s_int = "2\n2 2\n\n1 2 3 4 5 6 7 8"
            s_float = "2\n2 2\n\n1 2 3 4.5 5 6 7 8"
            s_exp = "2\n2 2\n\n1 2 3 4 5 6 7 1e3"

            @test parse_gam(s_int) isa NormalFormGame{2,Int}
            @test parse_gam(s_float) isa NormalFormGame{2,Float64}
            @test parse_gam(s_exp) isa NormalFormGame{2,Float64}
            @test parse_gam(s_exp)[2, 2][2] == 1000.0

            @test parse_gam(Float64, s_int) isa NormalFormGame{2,Float64}
            @test parse_gam(BigInt, s_int) isa NormalFormGame{2,BigInt}
            @test parse_gam(Rational{Int}, s_float) isa NormalFormGame{2,Rational{Int}}
            @test parse_gam(Rational{Int}, s_float)[2, 2][1] == 9//2
            @test_throws ArgumentError parse_gam(Int, s_float)

            g = NormalFormGame(Player([1 2; 3 4]), Player([5 6; 7 8]))
            @test gam_string(NormalFormGame(Float64, g)) ==
                  "2\n2 2\n\n1.0 3.0 2.0 4.0 5.0 6.0 7.0 8.0\n"

            # Integers that do not fit in Int
            s_big = "2\n2 2\n\n$(big(2)^70) 2 3 4 5 6 7 -$(big(2)^70)"
            @test parse_gam(s_big) isa NormalFormGame{2,BigInt}
            @test parse_gam(s_big)[1, 1][1] == big(2)^70
            @test gam_string(parse_gam(s_big)) == s_big * "\n"

            # Bool payoffs are written as 0 and 1
            g_bool = NormalFormGame(Player([true false; false true]),
                                    Player([true false; false true]))
            @test gam_string(g_bool) == "2\n2 2\n\n1 0 0 1 1 0 0 1\n"
            @test same_game(parse_gam(gam_string(g_bool)), g_bool)

            # Payoffs are not rounded when the caller's context is compact
            p_float = GAMPayoffVector((1, 1), [1.12341234, 2.0])
            @test sprint(write_gam, p_float; context=:compact => true) ==
                  "2\n1 1\n\n1.12341234 2.0\n"

            # An unsupported game is rejected before the file is opened
            p = GAMPayoffVector{2,UnsupportedReal}((1, 1), fill(UnsupportedReal(), 2))
            @test_throws MethodError write_gam(IOBuffer(), p)
            mktempdir() do dir
                path = joinpath(dir, "game.gam")
                write(path, "old content\n")
                @test_throws MethodError write_gam(path, p)
                @test read(path, String) == "old content\n"
            end
        end

        @testset "Whitespace" begin
            g = parse_gam("2\n3 2\n\n3 2 0 3 5 6 3 2 3 2 6 1")
            @test same_game(parse_gam("  2 3 2 3 2 0 3 5 6\n3 2 3 2 6 1\n\n"), g)
            @test same_game(parse_gam("2\r\n3 2\r\n\r\n3 2 0 3 5 6 3 2 3 2 6 1\r\n"), g)
        end

        @testset "File from QuantEcon.py" begin
            # Copied from quantecon/game_theory/tests/game_files in QuantEcon.py
            path = joinpath(@__DIR__, "game_files", "minimum_effort_game.gam")
            g = read_gam(path)

            @test g isa NormalFormGame{3,Float64}
            @test g.nums_actions == (3, 3, 3)
            @test g[1, 1, 1] == [1, 1, 1]
            @test g[1, 1, 3] == [1, 1, -19]
            @test g[2, 2, 2] == [2, 2, 2]
            @test g[3, 3, 3] == [3, 3, 3]
            @test same_game(parse_gam(gam_string(g)), g)
        end

        @testset "Invalid inputs" begin
            for s in ["", "  \n", "x", "0", "-1", "2\n3", string(typemax(Int)),
                      "2\n3 x\n\n1 2",
                      "2\n3 0\n\n", "2\n3 2\n\n1 2 3",
                      "2\n3 2\n\n1 2 3 4 5 6 7 8 9 10 11 12 13",
                      "2\n3 2\n\n1 2 3 4 5 6 7 8 9 10 11 z"]
                @test_throws ArgumentError parse_gam(s)
                @test_throws ArgumentError parse_gam(Float64, s)
            end
        end

        @testset "Number parsing and printing" begin
            _parse_exact = GameTheory._parse_exact
            _parse_payoffs = GameTheory._parse_payoffs

            for (tok, x) in [("3", 3), ("-3", -3), ("+3", 3), ("0.1", 1//10),
                             (".5", 1//2), ("-.5", -1//2), ("5.", 5),
                             ("-12.5e-3", -1//80),
                             ("+2e2", 200), ("1E-7", 1//10^7), ("1e30", big(10)^30),
                             ("1/3", 1//3), ("-1/3", -1//3), ("+1/3", 1//3),
                             ("6/4", 3//2)]
                @test _parse_exact(tok) == x
                @test _parse_exact(tok) isa Rational{BigInt}
            end
            for tok in ["", ".", "-", "1e", "1.2.3", "1/0", "1/2/3", "0.5/2",
                        "0x10", "abc", ".-5", ".+5", "1.-5"]
                @test_throws ArgumentError _parse_exact(tok)
            end

            @test _parse_payoffs(split("1 -2 +3")) isa Vector{Int}
            @test _parse_payoffs(split("1 1/3 -2")) isa Vector{Rational{BigInt}}
            @test _parse_payoffs(split("1 1/3 -2")) == [1, 1//3, -2]
            @test _parse_payoffs(split("1 0.5 -2")) isa Vector{Float64}
            @test _parse_payoffs(split("1 1e3 -2")) isa Vector{Float64}
            @test _parse_payoffs(split("1/3 0.5")) isa Vector{Float64}
            @test _parse_payoffs(split("1/3 0.5")) == [1/3, 0.5]

            # Integers that do not fit in Int
            x = big(10)^30 + 1
            @test _parse_payoffs(split("1 $x -3")) isa Vector{BigInt}
            @test _parse_payoffs(split("1 $x -3")) == [1, x, -3]
            for str in ["1 $x 1/3", "1/3 $x 1", "$x 1 1/3"]
                @test _parse_payoffs(split(str)) isa Vector{Rational{BigInt}}
                @test sort(_parse_payoffs(split(str))) == [1//3, 1, x]
            end
            @test _parse_payoffs(split("1 $x 0.5")) isa Vector{Float64}

            # What the writers print is read back without loss
            xs = [big(2)^70//1, 1//3, -5//2, big(3)^50//7]
            tokens = [sprint(GameTheory._print_payoff, x) for x in xs]
            @test _parse_payoffs(tokens) == xs
            @test _parse_payoffs(tokens) isa Vector{Rational{BigInt}}

            @test _parse_payoffs(Rational{Int}, split("0.1 1e-2 1/3")) ==
                  [1//10, 1//100, 1//3]
            @test _parse_payoffs(Float32, split("1/4 2")) == Float32[0.25, 2]
            @test_throws InexactError _parse_payoffs(Rational{Int}, ["1e30"])
            @test_throws ArgumentError _parse_payoffs(Int, ["1/3"])

            for (x, str) in [(1//3, "1/3"), (-5//2, "-5/2"), (3//1, "3"),
                             (big(10)^30//7, "1" * "0"^30 * "/7"), (-4, "-4"),
                             (0.25, "0.25"), (true, "1"), (false, "0")]
                @test sprint(GameTheory._print_payoff, x) == str
            end
        end

        @testset "Type inference" begin
            g = random_game(MersenneTwister(0), (3, 2))
            @inferred write_gam(IOBuffer(), g)
            @inferred gam_string(g)
            tokens = split("1 2 3 4.5")
            @inferred GameTheory._parse_payoffs(Float64, tokens)
            @inferred GameTheory._parse_payoffs(BigInt, split("1 2 3 4"))
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
