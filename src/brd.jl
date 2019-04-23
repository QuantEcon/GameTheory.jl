#=
    Tools for best response dynamics
=#

using StatsBase


# AbstractBRD

"""
    AbstractBRD

Abstract type representing the best response dynamics model.
"""
abstract type AbstractBRD{T<:Real} end

"""
    BRD

Type representing the best response dynamics model.

# Fields

- `N::Int` : The number of players.
- `player::Player{2,T}` : `Player` instance in the model.
- `num_actions::Int` : The number of actions for players.
"""
struct BRD{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
end

"""
    BRD(N, payoff_array)

Create a new BRD instance.

# Arguments

- `N::Integer` : The number of players.
- `payoff_array::Matrix` : Payoff array for each player.

# Returns

- `::BRD` : The best response dynamics model.
"""
function BRD(payoff_array::Matrix{T}, N::Integer) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return BRD(N, Player(payoff_array), num_actions)
end

"""
    KMR

Type representing the Kandori Mailath Rob model.

# Fields

- `N::Int` : The number of players.
- `player::Player` : `Player` instance in the model.
- `num_actions::Int` : The number of actions for players.
- `epsilon::Float64` : The probability of strategy flips.
"""
struct KMR{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
    epsilon::Float64
end

"""
    KMR(N, payoff_array, epsilon)

Create a new KMR instance.

# Arguments

- `N::Integer` : The number of players.
- `payoff_array::Matrix` : The payoff array for each player.
- `epsilon::Float64` : The probability of strategy flips.

# Returns

- `::KMR` : The Kandori Mailath Rob model.
"""
function KMR(payoff_array::Matrix{T},
             N::Integer,
             epsilon::Float64) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return KMR(N, Player(payoff_array), num_actions, epsilon)
end

"""
    SamplingBRD

Type representing the sampling best response dynamics model.

# Fields

- `N::Int` : The number of players.
- `player::Player` : `Player` instance in the model.
- `num_actions::Int` : The number of actions for players.
- `k::Int` : Sample size.
"""
struct SamplingBRD{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
    k::Int  #sample size
end

"""
    SamplingBRD(N, payoff_array, k)

Create a new SamplingBRD instance.

# Arguments

- `N::Integer` : The number of players.
- `payoff_array::Matrix` : Payoff array for a player.
- `k::Integer` : Sample size.

# Returns

- `::SamplingBRD` : The sampling best response dynamics model.
"""
function SamplingBRD(payoff_array::Matrix{T},
                     N::Integer,
                     k::Integer) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return SamplingBRD(N, Player(payoff_array), num_actions, k)
end


# play!

function play!(rng::AbstractRNG,
               brd::BRD,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions=BROptions())
    action_dist[action] -= 1
    next_action = best_response(brd.player, action_dist, options)
    action_dist[next_action] += 1
    return action_dist
end

function play!(rng::AbstractRNG,
               brd::KMR,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions=BROptions())
    action_dist[action] -= 1
    if rand(rng) <= brd.epsilon
        next_action = rand(rng, 1:brd.num_actions)
    else
        next_action = best_response(brd.player, action_dist, options)
    end
    action_dist[next_action] += 1
    return action_dist
end

function play!(rng::AbstractRNG,
               brd::SamplingBRD,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions=BROptions())
    action_dist[action] -= 1
    actions = sample(1:brd.num_actions, Weights(action_dist), brd.k)
    sample_action_dist = zeros(Int, brd.num_actions)
    for a in actions
        sample_action_dist[a] += 1
    end
    next_action = best_response(brd.player, sample_action_dist, options)
    action_dist[next_action] += 1
    return action_dist
end

@doc """
    play!([rng=Random.GLOBAL_RNG, ]brd, action, action_dist[, options=BROptions()])

Update an action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : `AbstractBRD` instance.
- `action::Integer` : A specified action.
- `action_dist::Vector{<:Integer}` : The distribution of players' actions.
- `options::BROptions` : Options for `best response` method.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""


# play

"""
    play([rng=Random.GLOBAL_RNG, ]brd, init_action_dist[, options=BROptions(); num_reps=1])

Return the action distribution after `num_reps` times iteration

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : `AbstractBRD` instance.
- `init_action_dist::Vector{<:Integer}` : The initial distribution of players' actions.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iterations.

# Returns

- `::Vector{<:Integer}` : The action distribution after iterations.
"""
function play(rng::AbstractRNG,
              brd::AbstractBRD,
              init_action_dist::Vector{<:Integer},
              options::BROptions=BROptions();
              num_reps::Integer=1)
    if length(init_action_dist) != brd.num_actions
        throw(ArgumentError("The length of init_action_dist must be the number
                             of actions"))
    end
    if sum(init_action_dist) != brd.N
        throw(ArgumentError("The sum of init_action_dist must be the number of
                             players"))
    end

    player_ind_seq = rand(rng, 1:brd.N, num_reps)
    for t in 1:num_reps
        action = searchsortedfirst(accumulate(+, init_action_dist),
                                   player_ind_seq[t])
        init_action_dist = play!(rng, brd, action, init_action_dist, options)
    end
    return init_action_dist
end

play(brd::AbstractBRD, init_action_dist::Vector{<:Integer},
     options::BROptions=BROptions(); num_reps::Integer=1) =
    play(Random.GLOBAL_RNG, brd, init_action_dist, options, num_reps=num_reps)


# time_series!

"""
    time_series!(rng, brd, out, player_ind_seq, options)

Update the matrix `out` which is used in `time_series` method given a player
index sequence.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : Instance of the model.
- `out::Matrix{<:Integer}` : Matrix representing the time series of action
  profiles.
- `player_ind_seq::Vector{<:Integer}` : The vector of player index.
- `options::BROptions` : Options for `best_response` method.

# Returns

- `out::Matrix{<:Integer}` : Updated `out`.
"""
function time_series!(rng::AbstractRNG,
                      brd::AbstractBRD,
                      out::Matrix{<:Integer},
                      player_ind_seq::Vector{<:Integer},
                      options::BROptions)
    ts_length = size(out, 2)
    action_dist = [out[i,1] for i in 1:brd.num_actions]
    for t in 1:ts_length-1
        action = searchsortedfirst(accumulate(+, action_dist), player_ind_seq[t])
        action_dist = play!(rng, brd, action, action_dist, options)
        for i in 1:brd.num_actions
            out[i,t+1] = action_dist[i]
        end
    end
    return out
end


# time_series

function time_series(rng::AbstractRNG,
                     brd::AbstractBRD,
                     ts_length::Integer,
                     init_action_dist::Vector{<:Integer},
                     options::BROptions=BROptions())
    if length(init_action_dist) != brd.num_actions
        throw(ArgumentError("The length of init_action_dist must be the number
                             of actions"))
    end
    if sum(init_action_dist) != brd.N
        throw(ArgumentError("The sum of init_action_dist must be the number of
                             players"))
    end
    
    player_ind_seq = rand(rng, 1:brd.N, ts_length)
    out = Matrix{Int}(undef, brd.num_actions, ts_length)
    for i in 1:brd.num_actions
        out[i, 1] = init_action_dist[i]
    end
    time_series!(rng, brd, out, player_ind_seq, options)
end

time_series(brd::AbstractBRD, ts_length::Integer,
            init_action_dist::Vector{<:Integer},
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, brd, ts_length, init_action_dist, options)

function time_series(rng::AbstractRNG,
                     brd::AbstractBRD,
                     ts_length::Integer,
                     options::BROptions=BROptions())
    player_ind_seq = rand(rng, 1:brd.N, ts_length)
    nums_actions = ntuple(i -> brd.num_actions, brd.N)
    init_actions = random_pure_actions(rng, nums_actions)
    action_dist = zeros(Int, brd.num_actions)
    for i in 1:brd.N
        action_dist[init_actions[i]] += 1
    end
    time_series(rng, brd, ts_length, action_dist, options)
end

time_series(brd::AbstractBRD, ts_length::Integer,
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, brd, ts_length, options)

@doc """
    time_series([rng=Random.GLOBAL_RNG, ]brd, ts_length, init_action_dist[, options=BROptions()])

Return the time series of action distribution.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : `AbstractBRD` instance.
- `ts_length::Integer` : The length of time series.
- `init_action_dist::Vector{<:Integer}` : Initial action distribution. If not
  provided, it is selected randomly.
- `options::BROptions` : Options for `best_response` method.

# Returns

- `::Matrix{<:Integer}` : The time series of action distributions.
"""
