#=
Benchmarks for generators/bimatrix_generators.jl

Each case times the construction of a game instance from the test suite
of von Stengel et al.; the random generators draw a fresh instance per
evaluation, advancing the case's own fixed-seed RNG.
=#
using GameTheory.Generators
using BenchmarkTools
using Random

#= Suite =#

suite = BenchmarkGroup()

# A fresh generator per case, so that adding or reordering cases does not
# alter the game data of the other cases
new_gen_rng() = MersenneTwister(0)

# blotto_game
suite["blotto_game"] = BenchmarkGroup()
rho = 0.5
for (h, t) in ((3, 62), (4, 21))
    rng = new_gen_rng()
    suite["blotto_game"]["h$(h)_t$(t)"] =
        @benchmarkable blotto_game($rng, $h, $t, $rho)
end

# ranking_game
let rng = new_gen_rng(), n = 2000
    suite["ranking_game"] = @benchmarkable ranking_game($rng, $n)
end

# sgc_game
let k = 500
    suite["sgc_game"] = @benchmarkable sgc_game($k)
end

# tournament_game
let rng = new_gen_rng(), n = 200, k = 2
    suite["tournament_game"] = @benchmarkable tournament_game($rng, $n, $k)
end

# unit_vector_game
suite["unit_vector_game"] = BenchmarkGroup()
for b in (true, false)
    rng = new_gen_rng()
    suite["unit_vector_game"]["avoid_pure_nash_$b"] =
        @benchmarkable unit_vector_game($rng, 2000; avoid_pure_nash=$b)
end

suite
