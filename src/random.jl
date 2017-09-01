#=
Generate random NormalFormGame instances.

Authors: Daisuke Oyama, Zejin Shi
=#

#
# Random Games Generating
#
"""
    random_game{N}([rng::AbstractRNG], nums_actions::NTuple{N,Int})

Return a random N-player NormalFormGame instance where the
payoffs are drawn independently from the uniform distribution
on [0, 1).

# Arguments

* `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
* `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions,
  one for each player.

# Returns

* `::NormalFormGame`: The generated random N-player NormalFormGame.
"""
function random_game{N}(rng::AbstractRNG, nums_actions::NTuple{N,Int})
    if N == 0
        throw(ArgumentError("nums_actions must be non-empty"))
    end

    players::NTuple{N,Player{N,Float64}} =
        ntuple(i -> Player(rand(rng, tuple(nums_actions[i:end]...,
                                      nums_actions[1:i-1]...))),
               N)

    return NormalFormGame(players)
end

random_game{N}(nums_actions::NTuple{N,Int}) = 
    random_game(Base.GLOBAL_RNG, nums_actions)

#
# Covariance Games Generating
#
"""
    covariance_game{N}([rng::AbstractRNG], nums_actions::NTuple{N,Int}, 
                       rho::Real)

Return a random N-player NormalFormGame instance with N>=2 where
the payoff profiles are drawn independently from the standard
multi-normal with the covariance of any pair of payoffs equal to
`rho`, as studied in Rinott and Scarsini (2000).

# Arguments

* `rng::AbstractRNG=GLOBAL_RNG`: Random number generator used.
* `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions, 
  one for each　player.
* `rho::T`: Covariance of a pair of payoff values. Must be in
  [-1/(N-1), 1], where N is the number of players.

# Returns

* `::NormalFormGame`: The generated random N-player NormalFormGame.

# References

* Y. Rinott and M. Scarsini, "On the Number of Pure Strategy
  Nash Equilibria in Random Games," Games and Economic Behavior
  (2000), 274-293.
"""
function covariance_game{N}(rng::AbstractRNG, nums_actions::NTuple{N,Int},
                            rho::Real)
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

    payoff_profile_array =
        reshape(transpose(x), (nums_actions..., N))

    return NormalFormGame(payoff_profile_array)
end

covariance_game{N}(nums_actions::NTuple{N,Int}, rho::Real) =
    covariance_game(Base.GLOBAL_RNG, nums_actions, rho)
