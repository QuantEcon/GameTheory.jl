using Games
using Games.Generators
using BenchmarkTools

const SUITE = BenchmarkGroup()

SUITE["support_enumeration"] = BenchmarkGroup(["support_enumeration"])
SUITE["support_enumeration"]["Float"] = BenchmarkGroup()
seed = 0
rng = MersenneTwister(seed)
ns = [10, 11]
for n in ns
    sz = (n, n)
    g = random_game(rng, sz)
    SUITE["support_enumeration"]["Float"][sz] =
        @benchmarkable support_enumeration($g)
end
SUITE["support_enumeration"]["Rational"] = BenchmarkGroup()
T = Rational{Int}
ns = [7, 8]
for n in ns
    sz = (n, n)
    g = NormalFormGame(eye(T, n))
    SUITE["support_enumeration"]["Rational"][sz] =
        @benchmarkable support_enumeration($g)
end


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
    @benchmarkable tournament_game($n, $k; seed=$seed)

# unit_vector_game
n = 2000
bools = [true, false]
SUITE["bimatrix_generators"]["unit_vector_game"] = BenchmarkGroup()
for b in bools
    SUITE["bimatrix_generators"]["unit_vector_game"][b] =
        @benchmarkable unit_vector_game($rng, $n; random=$b)
end
