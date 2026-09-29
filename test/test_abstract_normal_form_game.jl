# ------------------------------------ #
# Testing AbstractNormalFormGame        #
# ------------------------------------ #

using Distributions: Normal
using StaticArrays: SVector
using GameTheory: nums_actions, payoffs

# Games defined by their rules, at top level since types must be defined there

# 2 players, Int payoffs, `payoffs` returns a tuple
struct MatchingPenniesRules <: AbstractNormalFormGame{2,Int} end
GameTheory.nums_actions(::MatchingPenniesRules) = (2, 2)
GameTheory.payoffs(::MatchingPenniesRules, a) = a[1] == a[2] ? (1, -1) : (-1, 1)

# 3 players, Float64 payoffs, unequal numbers of actions, `payoffs` returns a
# vector
struct ThreePlayerRules <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::ThreePlayerRules) = (2, 3, 4)
GameTheory.payoffs(::ThreePlayerRules, a) =
    [a[1] + 10a[2], a[2] * a[3], (a[1] - a[3]) / 2]

# 1 player
struct OnePlayerRules <: AbstractNormalFormGame{1,Int} end
GameTheory.nums_actions(::OnePlayerRules) = (3,)
GameTheory.payoffs(::OnePlayerRules, a) = (a[1]^2,)

# `payoffs` of the wrong length
struct BadRules <: AbstractNormalFormGame{2,Int} end
GameTheory.nums_actions(::BadRules) = (2, 2)
GameTheory.payoffs(::BadRules, a) = (1, 2, 3)


@testset "Testing AbstractNormalFormGame" begin

    mp = MatchingPenniesRules()
    g_mp = NormalFormGame(Player([1 -1; -1 1]), Player([-1 1; 1 -1]))

    @testset "interface on NormalFormGame" begin
        @test NormalFormGame <: AbstractNormalFormGame
        @test @inferred(nums_actions(g_mp)) == (2, 2)
        @test @inferred(payoffs(g_mp, (1, 2))) == g_mp[1, 2]
        @test payoffs(g_mp, (1, 2)) isa SVector{2,Int}
        @test num_players(g_mp) == 2
    end

    @testset "SVector payoff profiles" begin
        @test @inferred(g_mp[1, 2]) isa SVector{2,Int}
        @test g_mp[1, 2] == [-1, 1]
        @test g_mp[1, 2] + g_mp[2, 1] == [-2, 2]
        @test payoff_profile_array(g_mp) isa Array{SVector{2,Int},2}
        g3 = NormalFormGame(ThreePlayerRules())
        @test g3[2, 3, 4] isa SVector{3,Float64}
        @test payoffs(OnePlayerRules(), (2,)) == (4,)
        @test payoffs(NormalFormGame(OnePlayerRules()), (2,)) isa SVector{1,Int}
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
            @test g3[Tuple(a)...] == payoffs(ThreePlayerRules(), Tuple(a))
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
        @test is_nash(OnePlayerRules(), 3)
        @test !is_nash(OnePlayerRules(), 1)
        @test is_nash(OnePlayerRules(), [0.0, 0.0, 1.0])
        @test payoff_profile_array(mp) == payoff_profile_array(g_mp)
        @test payoff_profile_array(tp) isa Array{SVector{3,Float64},3}
        @test payoff_profile_array(tp) == payoff_profile_array(g3)
        @test_throws MethodError delete_action(mp, 1, 2)
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

end
