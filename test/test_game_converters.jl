using GameTheory: PayoffProfileMatrix, PlayerMajor, ProfileMajor, as_vector
using LinearAlgebra: Transpose

using Random

# Element type that the writers do not accept
struct UnsupportedReal <: Real end

@testset "game_converters.jl" begin

    @testset "PayoffProfileMatrix" begin
        same_players(g1, g2) =
            all(g1.players[i].payoff_array == g2.players[i].payoff_array
                for i in 1:num_players(g1))

        @testset "Golden: N=3" begin
            nums_actions = (2, 3, 4)
            N = length(nums_actions)
            na = prod(nums_actions)

            A1 = reshape(collect(1:na), nums_actions)
            A2 = reshape(collect(101:100+na), nums_actions)
            A3 = reshape(collect(201:200+na), nums_actions)

            payoffs2d = hcat(vec(A1), vec(A2), vec(A3))
            payoffs_F = vec(payoffs2d)               # player-major
            payoffs_C = vec(permutedims(payoffs2d))  # profile-major

            payoffs4d = Array{Int,N+1}(undef, nums_actions..., N)
            payoffs4d[:, :, :, 1] .= A1
            payoffs4d[:, :, :, 2] .= A2
            payoffs4d[:, :, :, 3] .= A3
            g = NormalFormGame(payoffs4d)

            for (layout, payoffs1d) in [(PlayerMajor(), payoffs_F),
                                        (ProfileMajor(), payoffs_C)]
                p = @inferred PayoffProfileMatrix(nums_actions, payoffs1d, layout)

                @test p.nums_actions == nums_actions
                @test p.payoffs == payoffs2d
                @test num_players(p) == N
                @test as_vector(p, PlayerMajor()) == payoffs_F
                @test as_vector(p, ProfileMajor()) == payoffs_C

                g_from_p = @inferred NormalFormGame(p)
                @test g_from_p.nums_actions == g.nums_actions
                @test same_players(g_from_p, g)

                # Make an AbstractVector (SubArray) that equals payoffs1d
                payoffs1d_view = @view vcat([-999], payoffs1d, [999])[2:end-1]
                p = @inferred PayoffProfileMatrix(nums_actions, payoffs1d_view, layout)
                @test p.payoffs == payoffs2d
            end

            p = @inferred PayoffProfileMatrix(g)
            @test p.nums_actions == nums_actions
            @test p.payoffs == payoffs2d
            @test p.payoffs isa Matrix{Int}

            p = @inferred PayoffProfileMatrix(nums_actions, payoffs2d)
            @test p.payoffs === payoffs2d
            @test same_players(NormalFormGame(p), g)
        end

        @testset "Golden: N=2" begin
            # 3x2 game with payoff profiles, in column-major order over
            # (a_1, a_2):
            #   (1,1): (3,2)  (2,1): (0,6)  (3,1): (2,1)
            #   (1,2): (1,3)  (2,2): (4,0)  (3,2): (5,4)
            nums_actions = (3, 2)
            g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]))

            # Profile-major: (payoffs at profile 1)..., (payoffs at profile 2)...
            payoffs_C = [3, 2, 0, 6, 2, 1, 1, 3, 4, 0, 5, 4]
            # Player-major: (payoffs to player 1)..., (payoffs to player 2)...
            payoffs_F = [3, 0, 2, 1, 4, 5, 2, 6, 1, 3, 0, 4]

            p = PayoffProfileMatrix(g)
            @test as_vector(p, ProfileMajor()) == payoffs_C
            @test as_vector(p, PlayerMajor()) == payoffs_F

            for (layout, payoffs) in [(ProfileMajor(), payoffs_C),
                                      (PlayerMajor(), payoffs_F)]
                p = PayoffProfileMatrix(nums_actions, payoffs, layout)
                @test same_players(NormalFormGame(p), g)
            end
        end

        @testset "Round trip: N=$(length(ns))" for ns in [(4, 3), (2, 2, 3, 2)]
            N = length(ns)
            rng = MersenneTwister(12345)
            g = random_game(rng, 0:99, ns)
            p = @inferred PayoffProfileMatrix(g)
            g2 = @inferred NormalFormGame(p)

            p_BI = @inferred PayoffProfileMatrix(BigInt, g)
            g3 = @inferred NormalFormGame(Int, p_BI)

            for g_new in [g2, g3]
                @test g_new.nums_actions == g.nums_actions
                @test same_players(g_new, g)
            end

            # Through the flat vectors
            for layout in [PlayerMajor(), ProfileMajor()]
                v = @inferred as_vector(p, layout)
                @test v isa Vector{Int}
                p_new = PayoffProfileMatrix(ns, v, layout)
                @test p_new.payoffs == p.payoffs
                @test same_players(NormalFormGame(p_new), g)
            end
        end

        @testset "N=1" begin
            payoffs = [1., 2., 3.]
            nums_actions = (3,)

            g = NormalFormGame(Player(payoffs))
            ps = [PayoffProfileMatrix(nums_actions, payoffs, PlayerMajor()),
                  PayoffProfileMatrix(nums_actions, payoffs, ProfileMajor()),
                  PayoffProfileMatrix(g)]

            for p in ps
                @test p.nums_actions == nums_actions
                @test as_vector(p, PlayerMajor()) == payoffs
                @test as_vector(p, ProfileMajor()) == payoffs
                @test same_players(NormalFormGame(p), g)
            end
        end

        @testset "Storage and memory" begin
            nums_actions = (3, 2)
            v = collect(1:12)

            p_F = PayoffProfileMatrix(nums_actions, v, PlayerMajor())
            p_C = PayoffProfileMatrix(nums_actions, v, ProfileMajor())

            # Backed by the input vector, without copying
            @test p_F.payoffs isa Matrix{Int}
            @test p_C.payoffs isa Transpose{Int,Matrix{Int}}
            @test Base.mightalias(p_F.payoffs, v)
            @test Base.mightalias(p_C.payoffs, v)

            # as_vector shares memory if the layout matches the storage,
            # and copies otherwise
            for (p, matching, other) in [(p_F, PlayerMajor(), ProfileMajor()),
                                         (p_C, ProfileMajor(), PlayerMajor())]
                w = as_vector(p, matching)
                @test w isa Vector{Int}
                @test w == v
                @test Base.mightalias(w, v)

                w = as_vector(p, other)
                @test w isa Vector{Int}
                @test !Base.mightalias(w, v)

                # Player blocks are views; the game copies
                b = GameTheory._player_block(p, 2)
                @test size(b) == nums_actions
                @test Base.mightalias(b, p.payoffs)
                g = NormalFormGame(p)
                @test !Base.mightalias(g.players[2].payoff_array, p.payoffs)
            end

            # The player blocks read the same matrix in both storages
            for i in 1:2
                @test GameTheory._player_block(p_F, i) ==
                      reshape(p_F.payoffs[:, i], nums_actions)
                @test GameTheory._player_block(p_C, i) ==
                      reshape(p_C.payoffs[:, i], nums_actions)
            end
        end

        @testset "Invalid inputs" begin
            for layout in [PlayerMajor(), ProfileMajor()]
                @test_throws ArgumentError PayoffProfileMatrix((2, 2), [1, 2, 3], layout)
                @test_throws ArgumentError PayoffProfileMatrix((2, 0), Int[], layout)
            end
            @test_throws ArgumentError PayoffProfileMatrix((2, 2), zeros(Int, 2, 4))
            @test_throws ArgumentError PayoffProfileMatrix((2, 0), zeros(Int, 0, 2))
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

            v = [3, 0, 2, 1, 4, 5, 2, 6, 1, 3, 0, 4]
            p = PayoffProfileMatrix((3, 2), v, PlayerMajor())
            @test gam_string(p) == s
            # Written in player-major order whatever the storage
            p_C = PayoffProfileMatrix((3, 2), as_vector(p, ProfileMajor()),
                                      ProfileMajor())
            @test p_C.payoffs isa Transpose
            @test gam_string(p_C) == s
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
            p_float = PayoffProfileMatrix((1, 1), [1.12341234, 2.0], PlayerMajor())
            @test sprint(write_gam, p_float; context=:compact => true) ==
                  "2\n1 1\n\n1.12341234 2.0\n"

            # An unsupported game is rejected before the file is opened
            p = PayoffProfileMatrix((1, 1), fill(UnsupportedReal(), 2), PlayerMajor())
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

end
