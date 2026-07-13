#=
Benchmarks for lemke_howson.jl (the Lemke-Howson algorithm, with the
pivoting kernels of QuantEcon's pivoting.jl)

The plain cases run `lemke_howson` end to end on random games: n = 10 is
small enough that setup and allocations matter, while n = 100 is dominated
by the pivoting steps. Random games of size n = 200 can take hundreds of
thousands of pivoting steps from a fixed initial pivot, so that size is
run only with the capping heuristic of Codenotti et al. (`capping=10`),
which restarts with fresh initial pivots. The `_prealloc` case times the
repeated-solve regime through `lemke_howson!` with caller-owned arrays,
and `_full_workspace` additionally supplies `col_bufs` and `argmins` so
that the timed call performs no workspace allocations.
=#
using GameTheory
using BenchmarkTools
using Random

#= Suite =#

suite = BenchmarkGroup()

# A fresh generator per case, so that adding or reordering cases does not
# alter the game data of the other cases
new_lh_rng() = MersenneTwister(1234)

for n in (10, 100)
    g = random_game(new_lh_rng(), (n, n))
    # the fixture must converge to an equilibrium
    @assert lemke_howson(g, full_output=Val(true))[2].converged
    suite["random_n$n"] = @benchmarkable lemke_howson($g)
end

# Capping heuristic
for n in (100, 200)
    g = random_game(new_lh_rng(), (n, n))
    @assert lemke_howson(g, capping=10, full_output=Val(true))[2].converged
    suite["random_n$(n)_capping10"] =
        @benchmarkable lemke_howson($g, capping=10)
end

# Repeated-solve regime: caller-owned output and primary workspace arrays
# with default keywords (the auxiliary workspace is materialized lazily
# inside), and, in the second case, with the full workspace supplied so
# that the timed call performs no workspace allocations
let n = 10
    g = random_game(new_lh_rng(), (n, n))
    NE = (Vector{Float64}(undef, n), Vector{Float64}(undef, n))
    tableaux = (Matrix{Float64}(undef, n, 2n+1),
                Matrix{Float64}(undef, n, 2n+1))
    bases = (Vector{Int}(undef, n), Vector{Int}(undef, n))
    suite["random_n10_prealloc"] =
        @benchmarkable lemke_howson!($NE, $tableaux, $bases, $g)
    col_bufs = (Vector{Float64}(undef, n), Vector{Float64}(undef, n))
    argmins = Vector{Int}(undef, n)
    suite["random_n10_full_workspace"] =
        @benchmarkable lemke_howson!($NE, $tableaux, $bases, $g,
                                     col_bufs=$col_bufs, argmins=$argmins)
end

suite
