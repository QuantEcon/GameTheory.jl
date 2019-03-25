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

Type representing the best response dynamics model. Subtype of `AbstractBRD`.

# Fields

- `N::Int` : The number of players.
- `player::Player{2,T}` : Player instance in the model.
- `num_actions::Int` : The number of actions for each player.
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

- `::BRD`
"""
function BRD(N::Integer, payoff_array::Matrix{T}) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return BRD(N, Player(payoff_array), num_actions)
end

"""
    KMR

Type representing the Kandori Mailath Rob model. Subtype of `AbstractBRD`.

# Fields

- `N::Int` : The number of players.
- `player::Player` : Player instance in the model.
- `num_actions::Int` : The number of actions for a player.
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

- `::KMR`
"""
function KMR(N::Integer,
             payoff_array::Matrix{T},
             epsilon::Float64) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return BRD(N, Player(payoff_array), num_actions, epsilon)
end

"""
    SamplingBRD

Type representing sampling best response dynamics model.
Subtype of `AbstractBRD`.

# Fields

- `N::Int` : The number of players.
- `player::Player` : Player instance in the model.
- `num_actions::Int` : The number of actions for a player.
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

- `::SamplingBRD`
"""
function SamplingBRD(N::Integer,
                     payoff_array::Matrix{T},
                     k::Integer) where {T<:Real}
    num_actions = size(payoff_array, 1)
    if num_actions != size(payoff_array, 2)
        throw(ArgumentError("Payoff array must be square"))
    end
    return BRD(N, Player(payoff_array), num_actions, k)
end

# play!

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::BRD` : BRD instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best response` method.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
function play!(rng::AbstractRNG,
               brd::BRD,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions)
    action_dist[action] -= 1
    next_action = best_response(brd.player, action_dist, options)
    action_dist[next_action] += 1
    return action_dist
end

play!(brd::BRD, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::KMR` : KMR instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best response` method.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
function play!(rng::AbstractRNG,
               brd::KMR,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions)
    action_dist[action] -= 1
    if rand() <= brd.epsilon
        next_action = rand(rng, 1:brd.num_actions)
    else
        next_action = best_response(brd.player, actions, options)
    end
    action_dist[next_action] += 1
    return action_dist
end

play!(brd::KMR, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::SamplingBRD` : SamplingBRD instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best_response` method.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
function play!(rng::AbstractRNG,
               brd::SamplingBRD,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions)
    action_dist[action] -= 1
    actions = sample(1:brd.num_actions, Weights(action_dist), brd.k)
    sample_action_dist = zeros(brd.num_actions, dtype=Int)
    for a in actions
        sample_action_dist[a] += 1
    end
    next_action = best_response(brd.player, sample_action_dist, options)
    action_dist[next_action] += 1
    return action_dist
end

play!(brd::SamplingBRD, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)


# time_series!

"""
    time_series!(rng, brd, out, player_ind_seq)

Update `out`.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : Instance of the model.
- `out::Matrix{<:Integer}` : 
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
    ts_length = size(out, 1)
    action_dist = [out[i,1] for i in 1:brd.num_actions]
    for t in 1:ts_length
        action = searchsortedfirst(accumulate(+, action_dist), player_ind_seq[t])
        action_dist = play!(rng, brd, action, action_dist, options)
        for i in 1:brd.num_actions
            out[i,t+1] = action_dist[i]
        end
    end
    return out
end


# time_series

function _set_action_dist(brd::AbstractBRD, actions::PureActionProfile)
    if brd.N != length(actions)
        throw(ArgumentError("The length of action profile must
                             equal to the number of players"))
    end
    action_dist = zeros(brd.num_actions)
    for i in 1:brd.N
        action_dist[actions[i]] += 1
    end
    return action_dist
end

"""
    time_series(rng, brd, ts_length, init_actions, options)

Return the time series of action distribution.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : Instance of the model.
- `ts_length::Integer` : The length of time series.
- `init_actions::PureActionProfile` : Initial actions.
- `options::BROptions` : Options for `best_response` method;
    defaults to `BROptions()`.

# Returns

- `::Matrix{<:Integer}` : Time series of action distribution.
"""
function time_series(rng::AbstractRNG,
                     brd::AbstractBRD,
                     ts_length::Integer,
                     init_actions::PureActionProfile,
                     options::BROptions=BROptions())
    player_ind_seq = rand(1:brd.N, ts_length)
    action_dist = _set_action_dist(brd, init_actions)
    out = Matrix{Int}(undef, brd.num_actions, ts_length)
    for i in 1:brd.num_actions
        out[i,1] = action_dist[i]
    end
    time_series!(rng, brd, out, player_ind_seq, options)
end

time_series(brd::AbstractBRD, ts_length::Integer,
            init_actions::PureActionProfile, options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, brd, ts_length, init_actions, options)

"""
    time_series(rng, brd, ts_length, init_actions, options)

Return the time series of action distribution. Initial actions are choosed
randomly.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : Instance of the model.
- `ts_length::Integer` : The length of time series.
- `options::BROptions` : Options for `best_response` method;
    defaults to `BROptions()`.

# Returns

- `::Matrix{<:Integer}` : Time series of action distribution.
"""
function time_series(rng::AbstractRNG,
                     brd::AbstractBRD,
                     ts_length::Integer,
                     options::BROptions=BROptions())
    player_ind_seq = rand(1:brd.N, ts_length)
    nums_actions = ntuple(i -> brd.num_actions, brd.N)
    init_actions = random_pure_actions(nums_actions)
    time_series(rng, brd, ts_length, init_actions, options)
end

time_series(brd::AbstractBRD, ts_length::Integer,
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, brd, ts_length, options)
