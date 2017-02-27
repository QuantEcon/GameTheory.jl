#=
Generate random NormalFormGame instances.
=#

#
# Random Games Generating
#
"""
    random_game{N}(nums_actions::NTuple{N,Int})

Return a random NormalFormGame instance where the payoffs are drawn
independently from the uniform distribution on [0, 1).

# Arguements
* `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions,
    one for each player.

# Returns
* `::NormalFormGame`: The generated N-players random NormalFormGame.
"""
function random_game{N}(nums_actions::NTuple{N,Int})
    if N == 0
        throw(ArgumentError("nums_actions must be non-empty"))
    end

    players::NTuple{N,Player{N,Float64}} =
        ntuple(i -> Player(rand(tuple(nums_actions[i:end]...,
                                      nums_actions[1:i-1]...))),
               N)

    return NormalFormGame(players)
end

#
# Covariance Games Generating
#
"""
    covariance_game{N, T<:Real}(nums_actions::NTuple{N,Int}, rho::T)

Return a random NormalFormGame instance where the payoff profiles
are drawn independently from the standard multi-normal with the
covariance of any pair of payoffs equal to `rho`, as studied in
[1]_.

# Arguements
* `nums_actions::NTuple{N,Int}`: Tuple of the numbers of actions, 
    one for each　player.
* `rho::T`: Covariance of a pair of payoff values. Must be in
    [-1/(N-1), 1], where N is the number of players. T<:Real.

# Returns
* `::NormalFormGame`: The generated random NormalFormGame.

# References
1. Y. Rinott and M. Scarsini, "On the Number of Pure Strategy
   Nash Equilibria in Random Games," Games and Economic Behavior
   (2000), 274-293.
"""
function covariance_game{N, T<:Real}(nums_actions::NTuple{N,Int}, rho::T)
    if N <= 1
        throw(ArgumentError("length of nums_actions must be at least 2"))
    end

    if !(-1 / (N - 1) < rho < 1)
        lb = (N == 2) ? "-1" : "-1/$(N-1)"
        throw(ArgumentError("rho must be in ($lb, 1)"))
    end

    mu = zeros(N)
    C = fill(rho, (N, N))
    C[diagind(C)] = ones(N)

    d = MvNormal(mu, C)
    x = rand(d, prod(nums_actions))

    payoff_profile_array =
        reshape(transpose(x), (nums_actions..., N))

    return NormalFormGame(payoff_profile_array)
end
