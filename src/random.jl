#=
Generate random NormalFormGame instances.
=#

#
# Random Games Generating
#
"""
    random_game([rng=GLOBAL_RNG], [S=Float64], nums_actions)

Return a random N-player NormalFormGame instance where the payoffs are drawn
independently from the uniform distribution on the set as determined by `S`.
`S` is a range (such as `0:10`) or a subtype of `Integer` or `AbstractFloat`;
in the latter case, the set is [0, 1) for floats and `typemin(S):typemax(S)`
for integers.

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `S::Union{Type,AbstractRange}`: Set of values from which payoffs are drawn.
- `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions,
  one for each player.

# Returns

- `::NormalFormGame`: The generated random N-player NormalFormGame.

# Examples

```julia
julia> using GameTheory, Random

julia> rng = MersenneTwister(1234);

julia> g = random_game(rng, (4, 3));

julia> println(g)
4×3 NormalFormGame{2, Float64}:
 [0.590845, 0.066423]   [0.794026, 0.956753]  [0.246837, 0.646691]
 [0.766797, 0.112486]   [0.854147, 0.276021]  [0.579672, 0.651664]
 [0.566237, 0.0566425]  [0.200586, 0.842714]  [0.648882, 0.950498]
 [0.460085, 0.96467]    [0.298614, 0.945775]  [0.0109059, 0.789904]

julia> g = random_game(rng, 0:10, (4, 3));

julia> println(g)
4×3 NormalFormGame{2, Int64}:
 [6, 6]  [10, 10]  [5, 3]
 [8, 1]  [0, 6]    [3, 10]
 [9, 5]  [10, 8]   [4, 0]
 [9, 4]  [3, 1]    [9, 8]
```
"""
function random_game(rng::AbstractRNG, S::Union{Type{T},AbstractRange{T}},
                     nums_actions::NTuple{N,Int}) where {N,T<:Real}
    if N == 0
        throw(ArgumentError("nums_actions must be non-empty"))
    end

    players::NTuple{N,Player{N,T}} =
        ntuple(i -> Player(rand(rng, S, tuple(nums_actions[i:end]...,
                                              nums_actions[1:i-1]...))),
               N)

    return NormalFormGame(players)
end

random_game(rng::AbstractRNG, nums_actions::NTuple{N,Int}) where {N} =
    random_game(rng, Float64, nums_actions)

random_game(nums_actions::NTuple{N,Int}) where {N} =
    random_game(Random.GLOBAL_RNG, nums_actions)

random_game(S::Union{Type{T},AbstractRange{T}},
            nums_actions::NTuple{N,Int}) where {N,T<:Real} =
    random_game(Random.GLOBAL_RNG, S, nums_actions)

#
# Covariance Games Generating
#
"""
    covariance_game([rng=GLOBAL_RNG], nums_actions, rho)

Return a random N-player NormalFormGame instance with N>=2 where
the payoff profiles are drawn independently from the standard
multi-normal with the covariance of any pair of payoffs equal to
`rho`, as studied in Rinott and Scarsini (2000).

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions,
  one for each player.
- `rho::Real`: Covariance of a pair of payoff values. Must be in
  [-1/(N-1), 1], where N is the number of players.

# Returns

- `::NormalFormGame`: The generated random N-player NormalFormGame.

# References

- Y. Rinott and M. Scarsini, "On the Number of Pure Strategy
  Nash Equilibria in Random Games," Games and Economic Behavior
  (2000), 274-293.
"""
function covariance_game(rng::AbstractRNG, nums_actions::NTuple{N,Int},
                         rho::Real) where N
    if N <= 1
        throw(ArgumentError("length of nums_actions must be at least 2"))
    end

    if !(-1 / (N - 1) <= rho <= 1)
        lb = (N == 2) ? "-1" : "-1/$(N-1)"
        throw(ArgumentError("rho must be in [$lb, 1]"))
    end

    mu = zeros(N)
    Sigma = fill(rho, (N, N))
    Sigma[diagind(Sigma)] = ones(N)

    d = MVNSampler(mu, Sigma)
    x = rand(rng, d, prod(nums_actions))

    x_T = Matrix{eltype(x)}(undef, prod(nums_actions), N)
    transpose!(x_T, x)
    payoff_profile_array =
        reshape(x_T, (nums_actions..., N))

    return NormalFormGame(payoff_profile_array)
end

covariance_game(nums_actions::NTuple{N,Int}, rho::Real) where {N} =
    covariance_game(Random.GLOBAL_RNG, nums_actions, rho)

#
# Random action profile
#
"""
    random_pure_actions([rng=GLOBAL_RNG], nums_actions)

Return a tuple of random pure actions (integers).

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `nums_actions::NTuple{N,Int}`: N-tuple of the numbers of actions,
  one for each player.

# Returns

- `::NTuple{N,Int}`: N-tuple of random pure actions.

"""
random_pure_actions(rng::AbstractRNG, nums_actions::NTuple{N,Int}) where {N} =
    ntuple(i -> rand(rng, 1:nums_actions[i]), Val(N))

random_pure_actions(nums_actions::NTuple{N,Int}) where {N} =
    random_pure_actions(Random.GLOBAL_RNG, nums_actions)

"""
    random_mixed_actions([rng=GLOBAL_RNG], nums_actions)

Return a tuple of random mixed actions (vectors of floats).

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
- `nums_actions::NTuple{N,Int}`: N-tuple of the numbers of actions,
  one for each player.

# Returns

- `::NTuple{N,Vector{Float64}}`: N-tuple of random mixed actions.
"""
random_mixed_actions(rng::AbstractRNG, nums_actions::NTuple{N,Int}) where {N} =
    ntuple(i -> QuantEcon.random_probvec(rng, nums_actions[i]), Val(N))

random_mixed_actions(nums_actions::NTuple{N,Int}) where {N} =
    random_mixed_actions(Random.GLOBAL_RNG, nums_actions)
