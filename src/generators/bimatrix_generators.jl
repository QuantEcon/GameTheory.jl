#=
This module contains functions that generate NormalFormGame instances of the
2-player games studied by Fearnley, Igwe, and Savani (2015):

* Colonel Blotto Games (`blotto_game`): A non-zero sum extension of the Blotto
  game as studied by Hortala-Vallve and Llorente-Saguer (2012), where opposing
  parties have asymmetric and heterogeneous battlefield valuations.

* Ranking Games (`ranking_game`): These games were introduced by Goldberg,
  Goldberg, Krysta, and Ventre (2013) because of (i) they are
  important and well-motivated games by Economics literature and (ii) their
  computational complexity.

* SGC Games (`sgc_game`): These games were introduced by Sandholm, Gilpin, and
  Conitzer (2005) as a worst case scenario for support enumeration as it has a
  unique equilibrium where each player uses half of his actions in his support.

* Tournament Games (`tournament_game`): These games are constructed by
  Anbalagan et al. (2013) as games that do not have interim epsilon-Nash
  equilibria with constant cardinaliry supports for epsilon smaller than a
  certain threshold.

* Unit vector Games (`unit_vector_game`): These games were introduced by R. Savani
  and B. von Stengel (2015) whose payoffs for the column player are chosen randomly
  from the range [0, 1], but for the row player each column j contains exactly one
  1 payoff with the rest being 0.

Large part of the code here is based on the C code available at
https://github.com/bimatrix-games/bimatrix-generators distributed under BSD
3-Clause License.

References
----------
* Y. Anbalagan, S. Norin, R. Savani, and A. Vetta, "Polylogarithmic Supports
  Are Required for Approximate Well-Supported Nash Equilibria below 2/3," WINE,
  2013.

* J. Fearnley, T. P. Igwe, and R. Savani, "An Empirical Study of Finding
  Approximate Equilibria in Bimatrix Games," International Symposium on
  Experimental Algorithms (SEA), 2015.

* L.A. Goldberg, P.W. Goldberg, P. Krysta, and C. Ventre, "Ranking Games that
  have Competitiveness-based Strategies", Theoretical Computer Science, 2013

* R. Hortala-Vallve and A. Llorente-Saguer, "Pure Strategy Nash Equilibria in
  Non-Zero Sum Colonel Blotto Games", International Journal of Game Theory,
  2012.

* T. Sandholm, A. Gilpin, and V. Conitzer, "Mixed-Integer Programming Methods
  for Finding Nash Equilibria," AAAI, 2005.

* R. Savani and B. von Stengel, "Unit Vector Games," International
  Journal of Economic Theory, 2016.

=#

# BSD 3-Clause License
#
# Copyright (c) 2015, John Fearnley, Tobenna P. Igwe, Rahul Savani
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import QuantEcon: MVNSampler, simplex_grid, next_k_array!, k_array_rank
import LightGraphs: random_tournament_digraph


# blotto_game
"""
    blotto_game([rng=GLOBAL_RNG], h, t, rho[, mu=0])

Return a NormalFormGame instance of a 2-player non-zero sum Colonel Blotto game
(Hortala-Vallve and Llorente-Saguer, 2012), where the players have an equal
number `t` of troops to assign to `h` hills (so that the number of actions for
each player is equal to (t+h-1) choose (h-1) = (T+h-1)!/(T!*(h-1)!)). Each
player has a value for each hill that he receives if he assigns strictly more
troops to the hill than his opponent (ties are broken uniformly at random),
where the values are drawn from a multivariate normal distribution with
covariance `rho`. Each player’s payoff is the sum of the values of the hills won
by that player.


# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `h::Integer` : Number of hills.
- `t::Integer` : Number of troops.
- `rho::Real` : Covariance of the players' values of each hill. Must be in
  [-1, 1].
- `mu::Real=0` : Mean of the players' values of each hill.

# Returns

- `g::NormalFormGame`

# Examples

```julia
julia> rng = MersenneTwister(1234);

julia> g = blotto_game(rng, 2, 3, 0.5)
4×4 NormalFormGame{2,Float64}

julia> g.players[1]
4×4 Player{2,Float64}:
 0.186434  -0.494479  -0.494479  -0.494479
 0.867347   0.186434  -0.494479  -0.494479
 0.867347   0.867347   0.186434  -0.494479
 0.867347   0.867347   0.867347   0.186434

julia> g.players[2]
4×4 Player{2,Float64}:
 -0.688223  -1.02919   -1.02919   -1.02919
 -0.347259  -0.688223  -1.02919   -1.02919
 -0.347259  -0.347259  -0.688223  -1.02919
 -0.347259  -0.347259  -0.347259  -0.688223
```
"""
function blotto_game(rng::AbstractRNG, h::Integer, t::Integer, rho::Real,
                     mu::Real)
    d_mean = [mu; mu]
    d_cov = [1 rho; rho 1]
    dist = MVNSampler(d_mean, d_cov)
    values = rand(rng, dist, h)

    actions = simplex_grid(h, t)
    n = size(actions)[2]

    payoff_arrays = [Array{Float64}(n, n) for i in 1:2]

    payoffs = Array{Float64}(2)
    @inbounds for i = 1:n, j = 1:n
        fill!(payoffs, 0)
        for k = 1:h
            if actions[k, i] == actions[k, j]
                for p = 1:2
                    payoffs[p] += values[p, k] / 2
                end
            else
                winner = 1 + Int(actions[k, i] < actions[k, j])
                payoffs[winner] += values[winner, k]
            end
        end
        payoff_arrays[1][i, j] = payoffs[1]
        payoff_arrays[2][j, i] = payoffs[2]
    end

    g = NormalFormGame(
        [Player(payoff_array) for payoff_array in payoff_arrays]
    )

    return g
end

blotto_game(rng::AbstractRNG, h::Integer, t::Integer, rho::Real) =
    blotto_game(rng, h, t, rho, 0)

blotto_game(h::Integer, t::Integer, rho::Real, mu::Real) =
    blotto_game(Base.GLOBAL_RNG, h, t, rho, mu)

blotto_game(h::Integer, t::Integer, rho::Real) =
    blotto_game(Base.GLOBAL_RNG, h, t, rho, 0)

# ranking_game
"""
    ranking_game([rng=GLOBAL_RNG], k)

Return a NormalFormGame instance of the 2-player game introduced by Goldberg,
Goldberg, Krysta, and Ventre (2013) where each player chooses an effort level
associated with a score and a cost which are both increasing functions with
randomly generated step sizes.


# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `n::Integer` : Positive integer determining the number of actions, i.e,
   effort levels.

# Returns

- `g::NormalFormGame`

# Examples

```julia
julia> g = ranking_game(3)
3×3 NormalFormGame{2,Float64}

julia> g.players[1]
2×2 Player{2,Float64}:
 1.0       0.0        0.0
 0.9       0.9       -0.6
 0.633333  0.633333   0.633333

julia> g.players[2]
2×2 Player{2,Float64}:
 0.0        0.0         0.0
 0.7       -0.3        -0.3
 0.433333  -0.0666667  -0.566667
```
"""
function ranking_game(rng::AbstractRNG, n::Integer)
    score_p1, score_p2 = cumsum(rand(rng, 1:10, n)), cumsum(rand(rng, 1:10, n))

    cost_p1 = cumsum([0, rand(rng, 1:10, n-1)...])/(10*n)
    cost_p2 = cumsum([0, rand(rng, 1:10, n-1)...])/(10*n)

    tie_matrix = [s1==s2 for s1 in score_p1, s2 in score_p2]
    win_matrix_p1 = [s1>s2 for s1 in score_p1, s2 in score_p2]
    win_matrix_p2 = transpose(.!win_matrix_p1)
    payoff_arrays = [(win_matrix_p1+tie_matrix/2).-cost_p1,
                     (win_matrix_p2-transpose(tie_matrix)/2).-cost_p2]

    g = NormalFormGame(
        [Player(payoff_array) for payoff_array in payoff_arrays]
    )
    return g
end

ranking_game(n::Integer) = ranking_game(Base.GLOBAL_RNG, n)

# sgc_game
"""
    sgc_game(k)

Return a NormalFormGame instance of the 2-player game introduced by Sandholm,
Gilpin, and Conitzer (2005), which has a unique Nash equilibrium, where each
player plays half of the actions with positive probabilities. Payoffs are
normalized so that the minimum and the maximum payoffs are 0 and 1,
respectively.

# Arguments

- `k::Integer` : Positive integer determining the number of actions. The
  returned game will have `4*k-1` actions for each player.

# Returns

- `g::NormalFormGame`

# Examples

```julia
julia> g = sgc_game(2)
7×7 NormalFormGame{2,Float64}

julia> g.players[1]
7×7 Player{2,Float64}:
 0.75  0.5   1.0   0.5   0.5   0.5   0.5
 1.0   0.75  0.5   0.5   0.5   0.5   0.5
 0.5   1.0   0.75  0.5   0.5   0.5   0.5
 0.0   0.0   0.0   0.75  0.0   0.0   0.0
 0.0   0.0   0.0   0.0   0.75  0.0   0.0
 0.0   0.0   0.0   0.0   0.0   0.75  0.0
 0.0   0.0   0.0   0.0   0.0   0.0   0.75

julia> g.players[2]
7×7 Player{2,Float64}:
 0.75  0.5   1.0   0.5   0.5   0.5   0.5
 1.0   0.75  0.5   0.5   0.5   0.5   0.5
 0.5   1.0   0.75  0.5   0.5   0.5   0.5
 0.0   0.0   0.0   0.0   0.75  0.0   0.0
 0.0   0.0   0.0   0.75  0.0   0.0   0.0
 0.0   0.0   0.0   0.0   0.0   0.0   0.75
 0.0   0.0   0.0   0.0   0.0   0.75  0.0
```
"""
function sgc_game(k::Integer)
    n, m = 4*k-1, 2*k-1
    payoff_arrays = [Array{Float64}(n, n) for i in 1:2]

    for payoff_array in payoff_arrays
        for j in 1:m
            for i in 1:m
                payoff_array[i, j] = 0.75
            end
            for i in (m+1):n
                payoff_array[i, j] = 0
            end
        end
        for j in (m+1):n
            for i in 1:m
                payoff_array[i, j] = 0.5
            end
            for i in (m+1):n
                payoff_array[i, j] = 0
            end
        end

        payoff_array[1, m] = 1
        payoff_array[1, 2] = 0.5
        for i in 2:(m-1)
            payoff_array[i, i-1] = 1
            payoff_array[i, i+1] = 0.5
        end
        payoff_array[m, m-1] = 1
        payoff_array[m, 1] = 0.5
    end

    for h in 0:(k-1)
        i, j = m + 1 + 2*h, m + 1 + 2*h
        payoff_arrays[1][i, j] = payoff_arrays[1][i+1, j+1] = 0.75
        payoff_arrays[2][j, i+1] = payoff_arrays[2][j+1, j] = 0.75
    end

    g = NormalFormGame(
        [Player(payoff_array) for payoff_array in payoff_arrays]
    )
    return g
end

# tournament_game
"""
    tournament_game(n, k; seed=-1)

Return a NormalFormGame instance of the 2-player win-lose game, whose payoffs
are either 0 or 1, introduced by Anbalagan et al. (2013). Player 1 has n
actions, which constitute the set of nodes {1, ..., n}, while player 2 has n
choose k actions, each corresponding to a subset of k elements of the set of n
nodes. Given a randomly generated tournament graph on the n nodes, the payoff
for player 1 is 1 if, in the tournament, the node chosen by player 1 dominates
all the nodes in the k-subset chosen by player 2. The payoff for player 2 is 1
if player 2's k-subset contains player 1's node.

# Notes

The actions of player 2 are ordered according to the
[combinatorial number system]
(https://en.wikipedia.org/wiki/Combinatorial_number_system),
different from the order used in the original library in C.

# Arguments

- `n::Integer` : Number of nodes in the tournament graph.
- `k::Integer` : Size of subsets of nodes in the tournament graph.
- `seed::Integer=-1`: Seed for random number generator. If seed is negative,
  then `Base.GLOBAL_RNG` is used.

# Returns

- `g::NormalFormGame`

# Examples

```julia
julia> seed = 1234;

julia> g = tournament_game(5, 2; seed=seed)
5×10 NormalFormGame{2,Float64}

julia> g.players[1]
5×10 Player{2,Float64}:
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> g.players[2]
10×5 Player{2,Float64}:
 1.0  1.0  0.0  0.0  0.0
 1.0  0.0  1.0  0.0  0.0
 0.0  1.0  1.0  0.0  0.0
 1.0  0.0  0.0  1.0  0.0
 0.0  1.0  0.0  1.0  0.0
 0.0  0.0  1.0  1.0  0.0
 1.0  0.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0  1.0
 0.0  0.0  1.0  0.0  1.0
 0.0  0.0  0.0  1.0  1.0
```
"""
function tournament_game(n::Integer, k::Integer; seed::Integer=-1)

    m = zero(Csize_t)
    try
        m = binomial(Csize_t(n), Csize_t(k))
    catch InexactError
        error("Maximum allowed size exceeded")
    end

    R = zeros(Float64, n, m)
    C = zeros(Float64, m, n)
    tourn_graph = random_tournament_digraph(n; seed=seed)
    fadjlist = tourn_graph.fadjlist

    # populate matrix C
    X = collect(1:k)
    for j = 1:m
        C[j, X] = 1.
        next_k_array!(X)
    end

    # populate matrix R
    # continue to use array `X` to store indices
    a = Vector{Int}(k)
    for i = 1:n
        d = length(fadjlist[i])
        if d >= k
            for j = 1:k
                a[j] = j
            end
            while a[end] <= d
                for j in 1:k
                    X[j] = fadjlist[i][a[j]]
                end
                R[i, k_array_rank(Csize_t, X)] = 1.
                next_k_array!(a)
            end
        end
    end

    g = NormalFormGame([Player(R), Player(C)])

    return g
end

# unit_vector_game
"""
    unit_vector_game([rng=GLOBAL_RNG], k; random=true)

Return a NormalFormGame instance of the 2-player game introduced by R. Savani
and B. von Stengel (2015), which has the payoffs for the column player being
chosen randomly from the range [0, 1], and for the row player each column j
contains exactly one 1 payoff with the rest being 0.

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `k::Integer` : Positive integer determining the number of actions.
- `random::Bool=true` : If set `random` to be `false`, then 1 payoffs for
  the row player will be placed to avoid the existence of pure NE.

# Returns

- `g::NormalFormGame`

# Examples

```julia
julia> g = unit_vector_game(MersenneTwister(0), 5)
5×5 NormalFormGame{2,Float64}

julia> g.players[1]
5×5 Player{2,Float64}:
 0.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0

julia> g.players[2]
5×5 Player{2,Float64}:
 0.823648  0.203477   0.585812  0.655448  0.469304
 0.910357  0.0423017  0.539289  0.575887  0.0623676
 0.164566  0.0682693  0.260036  0.868279  0.353129
 0.177329  0.361828   0.910047  0.9678    0.767602
 0.27888   0.973216   0.167036  0.76769   0.043141

julia> pure_nash(g)
2-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (4, 4)

julia> g = unit_vector_game(MersenneTwister(0), 5, random=false)
5×5 NormalFormGame{2,Float64}

julia> g.players[1]
5×5 Player{2,Float64}:
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> g.players[2]
5×5 Player{2,Float64}:
 0.823648  0.203477   0.585812  0.655448  0.469304
 0.910357  0.0423017  0.539289  0.575887  0.0623676
 0.164566  0.0682693  0.260036  0.868279  0.353129
 0.177329  0.361828   0.910047  0.9678    0.767602
 0.27888   0.973216   0.167036  0.76769   0.043141

julia> pure_nash(g)
0-element Array{Tuple{Int64,Int64},1}
```
"""
function unit_vector_game(rng::AbstractRNG, k::Integer; random::Bool=true)

    payoff_arrays = [zeros(Float64, k, k) for i in 1:2]

    if random
        payoff_arrays[2][:, :] = rand(rng, (k, k))
        for c in 1:k
            payoff_arrays[1][rand(rng, 1:k), c] = 1
        end
    else
        while true
            payoff_arrays[2][:, :] = rand(rng, (k, k))
            brs = [indmax(payoff_arrays[2][:, r]) for r in 1:k]
            if length(unique(brs)) == 1
                continue
            end
            pure_nash_rows = [[] for i in 1:k]
            for (r, c) in enumerate(brs)
                push!(pure_nash_rows[c], r)
            end
            for c in 1:k
                payoff_arrays[1][rand(rng,
                                      deleteat!(collect(1:k),
                                                pure_nash_rows[c])),
                                 c] = 1
            end
            break
        end
    end

    g = NormalFormGame(
        [Player(payoff_array) for payoff_array in payoff_arrays]
    )

    return g

end

unit_vector_game(k::Integer; random::Bool=true) =
    unit_vector_game(Base.GLOBAL_RNG, k, random=random)
