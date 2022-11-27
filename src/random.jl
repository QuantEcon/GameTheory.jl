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
`S` is a range (such as `0:9`) or a subtype of `Integer` or `AbstractFloat`;
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

julia> rng = MersenneTwister(12345);

julia> g = random_game(rng, (4, 3));

julia> println(g)
4×3 NormalFormGame{2, Float64}:
 [0.562714, 0.586598]  [0.381128, 0.0501668]  [0.922317, 0.61179]
 [0.849939, 0.620099]  [0.365801, 0.215712]   [0.0404417, 0.569955]
 [0.371605, 0.965631]  [0.835014, 0.364706]   [0.573382, 0.923602]
 [0.283365, 0.754047]  [0.260024, 0.696476]   [0.981364, 0.0311643]

julia> g = random_game(rng, 0:9, (4, 3));

julia> println(g)
4×3 NormalFormGame{2, Int64}:
 [1, 5]  [1, 2]  [6, 2]
 [2, 5]  [0, 2]  [1, 0]
 [0, 5]  [3, 9]  [1, 1]
 [9, 5]  [2, 9]  [0, 6]
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

# Examples

```julia
julia> using GameTheory, Random

julia> rng = MersenneTwister(12345);

julia> g = covariance_game(rng, (4, 3), -0.7);

julia> println(g)
4×3 NormalFormGame{2, Float64}:
 [1.17236, -0.211696]   [1.46647, -1.13947]    [0.378353, 0.603951]
 [0.415565, 0.0779055]  [0.606808, 1.00812]    [1.12871, -1.03399]
 [0.685759, -0.278449]  [-0.588508, 0.464548]  [-0.970332, -0.0319236]
 [-1.47708, 1.12447]    [1.92585, -2.27959]    [-2.1476, 1.53569]
```

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
