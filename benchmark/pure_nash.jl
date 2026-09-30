#=
Benchmarks for pure_nash.jl (brute-force enumeration of pure-action Nash
equilibria)

`pure_nash` visits every action profile and checks each player's best
response, so its cost is the number of profiles times the cost of a
payoff-vector lookup per player. The cases are random games with 2, 3, and
4 players and roughly 10^4 profiles each; random games have few pure
equilibria, so the enumeration runs to the end.
=#
using GameTheory
using BenchmarkTools
using Random

#= Suite =#

suite = BenchmarkGroup()

# A fresh generator per case, so that adding or reordering cases does not
# alter the game data of the other cases
new_pn_rng() = MersenneTwister(0)

for nums_actions in ((100, 100), (20, 20, 20), (10, 10, 10, 10))
    g = random_game(new_pn_rng(), nums_actions)
    N = length(nums_actions)
    suite["random_$(N)p_n$(nums_actions[1])"] = @benchmarkable pure_nash($g)
end

suite
