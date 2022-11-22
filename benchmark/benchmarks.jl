using GameTheory
using BenchmarkTools
using Random
using LinearAlgebra

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
    g = NormalFormGame(Matrix{T}(I, n, n))
    SUITE["support_enumeration"]["Rational"][sz] =
        @benchmarkable support_enumeration($g)
end


SUITE["repeated_game"] = BenchmarkGroup(["repeated_game"])
pd_payoffs = [9.0 1.0; 10.0 3.0]
g = NormalFormGame(pd_payoffs)
rpd = RepeatedGame(g, 0.75)
SUITE["repeated_game"]["outerapproximation"] =
    @benchmarkable outerapproximation($rpd, nH=64, tol=1e-9)
