#=
Benchmarks for repeated_game.jl (outer approximation of the equilibrium
payoff set by the algorithm of Judd, Yeltekin, and Conklin)

The case runs `outerapproximation` on a prisoner's dilemma with discount
factor 0.75, with 64 subgradients; the iteration-to-iteration linear
programs dominate the cost.
=#
using GameTheory
using BenchmarkTools

#= Suite =#

suite = BenchmarkGroup()

pd_payoffs = [9.0 1.0; 10.0 3.0]
g = NormalFormGame(pd_payoffs)
rpd = RepeatedGame(g, 0.75)
suite["outerapproximation_nH64"] =
    @benchmarkable outerapproximation($rpd, nH=64, tol=1e-9)

suite
