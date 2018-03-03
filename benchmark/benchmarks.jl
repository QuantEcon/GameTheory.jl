using Games
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
