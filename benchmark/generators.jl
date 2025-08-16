using GameTheory.Generators
using BenchmarkTools
using Random

const SUITE = BenchmarkGroup()

SUITE["bimatrix_generators"] = BenchmarkGroup(["bimatrix_generators"])
seed = 0
rng = MersenneTwister(seed)

# blotto_game
hts = [(3, 62), (4, 21)]  # (h, t)
rho = 0.5
SUITE["bimatrix_generators"]["blotto_game"] = BenchmarkGroup()
for ht in hts
    SUITE["bimatrix_generators"]["blotto_game"][ht] =
        @benchmarkable blotto_game($rng, ($ht)..., $rho)
end

# ranking_game
n = 2000
SUITE["bimatrix_generators"]["ranking_game"] =
    @benchmarkable ranking_game($rng, $n)

# sgc_game
k = 500
SUITE["bimatrix_generators"]["sgc_game"] = @benchmarkable sgc_game($k)

# tournament_game
n, k = 200, 2
SUITE["bimatrix_generators"]["tournament_game"] =
    @benchmarkable tournament_game($rng, $n, $k)

# unit_vector_game
n = 2000
bools = [true, false]
SUITE["bimatrix_generators"]["unit_vector_game"] = BenchmarkGroup()
for b in bools
    SUITE["bimatrix_generators"]["unit_vector_game"][b] =
        @benchmarkable unit_vector_game($rng, $n; avoid_pure_nash=$b)
end
