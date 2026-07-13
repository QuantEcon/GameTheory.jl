#=
Benchmarks for support_enumeration.jl

The `Float` cases run on random games with `Float64` payoffs (LU-based
linear solves); the `Rational` cases run on identity-matrix games with
`Rational{Int}` payoffs, exercising the exact-arithmetic path. Sizes are
kept small since the cost grows combinatorially in the number of actions.
=#
using GameTheory
using BenchmarkTools
using Random
using LinearAlgebra

#= Suite =#

suite = BenchmarkGroup()

# A fresh generator per case, so that adding or reordering cases does not
# alter the game data of the other cases
new_se_rng() = MersenneTwister(0)

suite["Float"] = BenchmarkGroup()
for n in (10, 11)
    g = random_game(new_se_rng(), (n, n))
    suite["Float"]["random_n$n"] = @benchmarkable support_enumeration($g)
end

suite["Rational"] = BenchmarkGroup()
for n in (7, 8)
    g = NormalFormGame(Matrix{Rational{Int}}(I, n, n))
    suite["Rational"]["identity_n$n"] = @benchmarkable support_enumeration($g)
end

suite
