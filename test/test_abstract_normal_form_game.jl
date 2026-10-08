# ------------------------------------ #
# Testing AbstractNormalFormGame        #
# ------------------------------------ #

using Distributions: Normal
using StaticArrays: SVector
using GameTheory: nums_actions, payoff_profile, _payoff_profile, _payoff_vector

# Games defined by their rules, at top level since types must be defined there

# 2 players, Int payoffs, `payoff_profile` returns a tuple
struct MatchingPenniesRules <: AbstractNormalFormGame{2,Int} end
GameTheory.nums_actions(::MatchingPenniesRules) = (2, 2)
GameTheory.payoff_profile(::MatchingPenniesRules, a) =
    a[1] == a[2] ? (1, -1) : (-1, 1)

# 3 players, Float64 payoffs, unequal numbers of actions, `payoff_profile`
# returns a vector
struct ThreePlayerRules <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::ThreePlayerRules) = (2, 3, 4)
GameTheory.payoff_profile(::ThreePlayerRules, a) =
    [a[1] + 10a[2], a[2] * a[3], (a[1] - a[3]) / 2]

# 3 players, Float64 payoffs, `payoff_profile` returns a tuple of mixed element
# types; the same payoffs as `ThreePlayerRules`
struct MixedTupleRules <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::MixedTupleRules) = (2, 3, 4)
GameTheory.payoff_profile(::MixedTupleRules, a) =
    (a[1] + 10a[2], a[2] * a[3], (a[1] - a[3]) / 2)

# Matching pennies counting the calls to `payoff_profile`
struct CountingRules <: AbstractNormalFormGame{2,Int}
    calls::Base.RefValue{Int}
end
CountingRules() = CountingRules(Ref(0))
GameTheory.nums_actions(::CountingRules) = (2, 2)
function GameTheory.payoff_profile(g::CountingRules, a)
    g.calls[] += 1
    return a[1] == a[2] ? (1, -1) : (-1, 1)
end

# 1 player
struct OnePlayerRules <: AbstractNormalFormGame{1,Int} end
GameTheory.nums_actions(::OnePlayerRules) = (3,)
GameTheory.payoff_profile(::OnePlayerRules, a) = (a[1]^2,)

# `payoff_profile` of the wrong length
struct BadRules <: AbstractNormalFormGame{2,Int} end
GameTheory.nums_actions(::BadRules) = (2, 2)
GameTheory.payoff_profile(::BadRules, a) = (1, 2, 3)


@testset "Testing AbstractNormalFormGame" begin

    mp = MatchingPenniesRules()
    g_mp = NormalFormGame(Player([1 -1; -1 1]), Player([-1 1; 1 -1]))

    @testset "interface on NormalFormGame" begin
        @test NormalFormGame <: AbstractNormalFormGame
        @test @inferred(nums_actions(g_mp)) == (2, 2)
        @test @inferred(payoff_profile(g_mp, (1, 2))) == g_mp[1, 2]
        @test payoff_profile(g_mp, (1, 2)) isa SVector{2,Int}
        @test num_players(g_mp) == 2
    end

    @testset "SVector payoff profiles" begin
        @test @inferred(g_mp[1, 2]) isa SVector{2,Int}
        @test g_mp[1, 2] == [-1, 1]
        @test g_mp[1, 2] + g_mp[2, 1] == [-2, 2]
        @test payoff_profile_array(g_mp) isa Array{SVector{2,Int},2}
        g3 = NormalFormGame(ThreePlayerRules())
        @test g3[2, 3, 4] isa SVector{3,Float64}
        @test payoff_profile(OnePlayerRules(), (2,)) == (4,)
        @test payoff_profile(NormalFormGame(OnePlayerRules()), (2,)) isa
            SVector{1,Int}
    end

    @testset "payoff profile normalization" begin
        @test @inferred(_payoff_profile(mp, (1, 2))) === SVector(-1, 1)
        @test _payoff_profile(ThreePlayerRules(), (1, 1, 1)) isa
            SVector{3,Float64}
        @test_throws DimensionMismatch _payoff_profile(BadRules(), (1, 1))
        @test_throws DimensionMismatch payoff_profile_array(BadRules())
        @test_throws DimensionMismatch is_nash(BadRules(), (1, 1))

        mt = MixedTupleRules()
        tp = ThreePlayerRules()
        @test payoff_profile(mt, (1, 2, 3)) isa Tuple{Int,Int,Float64}
        @test @inferred(_payoff_profile(mt, (1, 2, 3))) isa SVector{3,Float64}
        g_mt = @inferred NormalFormGame(mt)
        g_tp = NormalFormGame(tp)
        @test g_mt isa NormalFormGame{3,Float64}
        for i in 1:3
            @test g_mt.players[i].payoff_array == g_tp.players[i].payoff_array
        end
        @test @inferred(payoff_profile_array(mt)) == payoff_profile_array(tp)
        for a in CartesianIndices(nums_actions(mt))
            @test is_nash(mt, Tuple(a)) == is_nash(g_mt, Tuple(a))
            @test is_pareto_efficient(mt, Tuple(a)) ==
                is_pareto_efficient(g_mt, Tuple(a))
            @test is_pareto_dominant(mt, Tuple(a)) ==
                is_pareto_dominant(g_mt, Tuple(a))
        end
        x = ([0.5, 0.5], [0.2, 0.3, 0.5], [0.25, 0.25, 0.25, 0.25])
        @test is_nash(mt, x) == is_nash(g_mt, x)
    end

    @testset "indexing and summary" begin
        @test @inferred(mp[1, 2]) === SVector(-1, 1)
        @test mp[1, 2] == g_mp[1, 2]
        @test mp[CartesianIndex(2, 1)] == g_mp[CartesianIndex(2, 1)]
        tp = ThreePlayerRules()
        g3 = NormalFormGame(tp)
        for a in CartesianIndices(nums_actions(tp))
            @test tp[Tuple(a)...] == g3[a]
        end
        @test @inferred(tp[2, 3, 4]) isa SVector{3,Float64}
        @test @inferred(OnePlayerRules()[2]) === 4
        @test NormalFormGame(OnePlayerRules())[2] === 4
        @test_throws DimensionMismatch mp[1]
        @test_throws DimensionMismatch mp[1, 2, 3]
        @test_throws DimensionMismatch tp[1, 2]
        @test summary(mp) == "2×2 MatchingPenniesRules"
        @test summary(tp) == "2×3×4 ThreePlayerRules"
        @test summary(OnePlayerRules()) == "3-element OnePlayerRules"
        @test summary(g_mp) == "2×2 NormalFormGame{2, Int64}"
    end

    @testset "tabulation" begin
        @test num_players(mp) == 2
        g = @inferred NormalFormGame(mp)
        @test g isa NormalFormGame{2,Int}
        for i in 1:2
            @test g.players[i].payoff_array == g_mp.players[i].payoff_array
        end

        g3 = @inferred NormalFormGame(ThreePlayerRules())
        @test g3 isa NormalFormGame{3,Float64}
        @test nums_actions(g3) == (2, 3, 4)
        for a in CartesianIndices(nums_actions(g3))
            @test g3[Tuple(a)...] ==
                payoff_profile(ThreePlayerRules(), Tuple(a))
        end

        g1 = @inferred NormalFormGame(OnePlayerRules())
        @test g1 isa NormalFormGame{1,Int}
        @test [g1[a] for a in 1:3] == [1, 4, 9]

        @test_throws DimensionMismatch NormalFormGame(BadRules())
    end

    @testset "eltype" begin
        g_rat = NormalFormGame(Rational{Int}, mp)
        @test g_rat isa NormalFormGame{2,Rational{Int}}
        @test g_rat[1, 1] == [1//1, -1//1]
        @test NormalFormGame(Float64, ThreePlayerRules()) isa
            NormalFormGame{3,Float64}
    end

    @testset "convert" begin
        @test convert(NormalFormGame, g_mp) === g_mp
        @test convert(NormalFormGame{2,Int}, g_mp) === g_mp
        @test NormalFormGame(g_mp) !== g_mp
        g = convert(NormalFormGame, mp)
        @test g isa NormalFormGame{2,Int}
        @test convert(NormalFormGame{2,Float64}, mp) isa NormalFormGame{2,Float64}
    end

    @testset "game functions on the interface" begin
        tp = ThreePlayerRules()
        g3 = NormalFormGame(tp)
        for a in CartesianIndices(nums_actions(g3))
            @test is_nash(tp, Tuple(a)) == is_nash(g3, Tuple(a))
            @test is_pareto_efficient(tp, Tuple(a)) ==
                is_pareto_efficient(g3, Tuple(a))
            @test is_pareto_dominant(tp, Tuple(a)) ==
                is_pareto_dominant(g3, Tuple(a))
        end
        @test is_nash(mp, ([0.5, 0.5], [0.5, 0.5]))
        @test !is_nash(mp, ([0.6, 0.4], [0.5, 0.5]))
        x = ([0.5, 0.5], [0.2, 0.3, 0.5], [0.25, 0.25, 0.25, 0.25])
        @test is_nash(tp, x) == is_nash(g3, x)
        @test_throws ArgumentError is_nash(tp, (1, 2))
        @test_throws ArgumentError is_nash(tp, (1, 2, 3, 4))
        @test_throws ArgumentError is_nash(tp, ([0.5, 0.5], [0.2, 0.3, 0.5]))
        @test_throws ArgumentError is_nash(tp, ())
        @test is_nash(OnePlayerRules(), 3)
        @test !is_nash(OnePlayerRules(), 1)
        @test is_nash(OnePlayerRules(), [0.0, 0.0, 1.0])
        @test payoff_profile_array(mp) == payoff_profile_array(g_mp)
        @test payoff_profile_array(tp) isa Array{SVector{3,Float64},3}
        @test payoff_profile_array(tp) == payoff_profile_array(g3)
        @test_throws MethodError delete_action(mp, 1, 2)
    end

    @testset "_payoff_vector" begin
        tp = ThreePlayerRules()
        g3 = NormalFormGame(tp)
        xs = ([0.5, 0.5], [0.2, 0.3, 0.5], [0.25, 0.25, 0.25, 0.25])
        for i in 1:3
            player = g3.players[i]
            for a in CartesianIndices(nums_actions(tp))
                opp = GameTheory.get_opponents_actions(Tuple(a), i)
                @test @inferred(_payoff_vector(tp, i, opp)) ==
                    payoff_vector(player, opp)
                @test _payoff_vector(g3, i, opp) == payoff_vector(player, opp)
            end
            opp = GameTheory.get_opponents_actions(xs, i)
            @test _payoff_vector(tp, i, opp) ≈ payoff_vector(player, opp)
            @test _payoff_vector(g3, i, opp) == payoff_vector(player, opp)
        end
        x = [0.3, 0.7]
        for i in 1:2
            player = g_mp.players[i]
            for a in 1:2
                @test _payoff_vector(mp, i, (a,)) == payoff_vector(player, a)
                @test _payoff_vector(g_mp, i, (a,)) == payoff_vector(player, a)
            end
            @test _payoff_vector(mp, i, (x,)) ≈ payoff_vector(player, x)
            @test _payoff_vector(g_mp, i, (x,)) == payoff_vector(player, x)
        end
        one = OnePlayerRules()
        @test _payoff_vector(one, 1, ()) == [1, 4, 9]
        @test _payoff_vector(NormalFormGame(one), 1, ()) == [1, 4, 9]
    end

    @testset "no tabulation in the functions on the interface" begin
        g = CountingRules()
        @test !is_nash(g, (1, 2))
        @test g.calls[] == 2 + 2             # n_i calls per player
        g = CountingRules()
        @test is_nash(g, ([0.5, 0.5], [0.5, 0.5]))
        @test g.calls[] == 2 * 2             # one pass over the profiles
        g = CountingRules()
        @test _payoff_vector(g, 1, (2,)) == [-1, 1]
        @test g.calls[] == 2
        g = CountingRules()
        @test g[2, 1] == [-1, 1]
        @test g.calls[] == 1
    end

    @testset "PayoffVector from the interface" begin
        tp = ThreePlayerRules()
        g3 = NormalFormGame(tp)
        for PV in (GameTheory.GAMPayoffVector, GameTheory.NFGPayoffVector)
            p = @inferred PV(tp)
            @test p.nums_actions == (2, 3, 4)
            @test p.payoffs == PV(g3).payoffs
            @test PV(MixedTupleRules()).payoffs == PV(g3).payoffs
            @test PV(Float32, tp).payoffs == PV(Float32, g3).payoffs
            @test_throws DimensionMismatch PV(BadRules())
        end
    end

    @testset "Nash equilibrium solvers" begin
        @test pure_nash(mp) == pure_nash(g_mp)
        @test pure_nash(ThreePlayerRules()) ==
            pure_nash(NormalFormGame(ThreePlayerRules()))
        @test lemke_howson(mp) == lemke_howson(g_mp)
        @test support_enumeration(mp) == support_enumeration(g_mp)
        @test vertex_enumeration(mp) == vertex_enumeration(g_mp)
        @test lrsnash(mp) == lrsnash(g_mp)
        NEs = hc_solve(mp, show_progress=false)
        @test length(NEs) == 1
        @test is_nash(g_mp, NEs[1])
    end

    @testset "learning algorithms, repeated game, converters" begin
        @test FictitiousPlay(mp).players[2].payoff_array ==
            FictitiousPlay(g_mp).players[2].payoff_array
        @test StochasticFictitiousPlay(mp, Normal()).players[1].payoff_array ==
            g_mp.players[1].payoff_array
        adj = [0 1; 1 0]
        @test LocalInteraction(mp, adj).players[1].payoff_array ==
            LocalInteraction(g_mp, adj).players[1].payoff_array
        @test LogitDynamics(mp, 1.0).players[1].payoff_array ==
            LogitDynamics(g_mp, 1.0).players[1].payoff_array

        @test RepeatedGame(mp, 0.5).sg.players[1].payoff_array ==
            g_mp.players[1].payoff_array

        @test GameTheory.GAMPayoffVector(mp).payoffs ==
            GameTheory.GAMPayoffVector(g_mp).payoffs
        @test GameTheory.GAMPayoffVector(Float64, mp).payoffs ==
            GameTheory.GAMPayoffVector(Float64, g_mp).payoffs
        @test gam_string(mp) == gam_string(g_mp)
        @test sprint(write_gam, mp) == gam_string(g_mp)
        mktempdir() do dir
            path = joinpath(dir, "mp.gam")
            write_gam(path, mp)
            @test read(path, String) == gam_string(g_mp)
        end

        @test GameTheory.NFGPayoffVector(mp).payoffs ==
            GameTheory.NFGPayoffVector(g_mp).payoffs
        @test nfg_string(mp) == nfg_string(g_mp)
        @test sprint(write_nfg, mp) == nfg_string(g_mp)
        mktempdir() do dir
            path = joinpath(dir, "mp.nfg")
            write_nfg(path, mp)
            @test read(path, String) == nfg_string(g_mp)
        end
    end

end
