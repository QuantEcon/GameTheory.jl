#=
This module contains functions that generate NormalFormGame instances of the
2-player games studied by Fearnley, Igwe, and Savani (2015):

* Colonel Blotto Games (`blotto_game`)

* Ranking Games (`ranking_game`)

* SGC Games (`sgc_game`): These games were introduced by Sandholm, Gilpin, and
  Conitzer (2005) as a worst case scenario for support enumeration as it has a
  unique equilibrium where each player uses half of his actions in his support.

* Tournament Games (`tournament_game`)

* Unit vector Games (`unit_vector_game`): These games were introduced by R. Savani
  and B. von Stengel (2015) whose payoffs for the column player are chosen randomly
  from the range [0, 1], but for the row player each column j contains exactly one
  1 payoff with the rest being 0.

Large part of the code here is based on the C code available at
https://github.com/bimatrix-games/bimatrix-generators distributed under BSD
3-Clause License.

References
----------
* J. Fearnley, T. P. Igwe, R. Savani, "An Empirical Study of Finding
  Approximate Equilibria in Bimatrix Games," International Symposium on
  Experimental Algorithms (SEA), 2015.

* T. Sandholm, A. Gilpin, and V. Conitzer, "Mixed-Integer Programming Methods
  for Finding Nash Equilibria," AAAI, 2005.

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
