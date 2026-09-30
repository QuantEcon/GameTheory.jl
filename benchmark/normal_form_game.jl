#=
Benchmarks for normal_form_game.jl (the `AbstractNormalFormGame` interface)

The generic code for a rules-defined game calls `payoffs` once per action
profile and normalizes the returned value to an `SVector` in
`_payoff_profile`. The cases are 3-player games with 30 actions each, about
2.7 * 10^4 profiles, whose `payoffs` return a tuple of `Float64`s, a tuple
of mixed element types, or a `Vector`. The operations are the tabulation
into a `NormalFormGame`, `payoff_profile_array`, and `is_nash` with a mixed
action profile, each a single pass over the profiles.
=#
using GameTheory
using BenchmarkTools

#= Games =#

const NUMS_ACTIONS_ANFG = (30, 30, 30)

struct HomTupleGame <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::HomTupleGame) = NUMS_ACTIONS_ANFG
GameTheory.payoffs(::HomTupleGame, a) =
    (a[1] + 0.5a[2], float(a[2] - a[3]), 0.1a[3] * a[1])

struct MixedTupleGame <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::MixedTupleGame) = NUMS_ACTIONS_ANFG
GameTheory.payoffs(::MixedTupleGame, a) =
    (a[1] + 0.5a[2], a[2] - a[3], 0.1a[3] * a[1])

struct VectorGame <: AbstractNormalFormGame{3,Float64} end
GameTheory.nums_actions(::VectorGame) = NUMS_ACTIONS_ANFG
GameTheory.payoffs(::VectorGame, a) =
    [a[1] + 0.5a[2], a[2] - a[3], 0.1a[3] * a[1]]

#= Suite =#

suite = BenchmarkGroup()

for op in ("NormalFormGame", "payoff_profile_array", "is_nash_mixed")
    suite[op] = BenchmarkGroup()
end

let x = ntuple(_ -> fill(1/30, 30), 3)
    for (name, g) in (("hom_tuple", HomTupleGame()),
                      ("mixed_tuple", MixedTupleGame()),
                      ("vector", VectorGame()))
        suite["NormalFormGame"][name] = @benchmarkable NormalFormGame($g)
        suite["payoff_profile_array"][name] =
            @benchmarkable payoff_profile_array($g)
        suite["is_nash_mixed"][name] = @benchmarkable is_nash($g, $x)
    end
end

suite
