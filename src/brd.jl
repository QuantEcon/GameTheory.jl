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
function BRD(payoff_array::Matrix{T}, N::Integer) where {T<:Real}
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

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::BRD` : BRD instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best response` method;
    defaults to `BROptions()`.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
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

play!(brd::BRD, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions=BROptions()) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::KMR` : KMR instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best response` method;
    defaults to `BROptions()`.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
function play!(rng::AbstractRNG,
               brd::KMR,
               action::Integer,
               action_dist::Vector{<:Integer},
               options::BROptions=BROptions())
    action_dist[action] -= 1
    if rand() <= brd.epsilon
        next_action = rand(rng, 1:brd.num_actions)
    else
        next_action = best_response(brd.player, action_dist, options)
    end
    action_dist[next_action] += 1
    return action_dist
end

play!(brd::KMR, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions=BROptions()) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)

"""
    play!(rng, brd, action, action_dist, options)

Update action distribution given a specified action.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::SamplingBRD` : SamplingBRD instance.
- `action::Integer` : Played action.
- `action_dist::Vector{<:Integer}` : The distribution of actions of players.
- `options::BROptions` : Options for `best_response` method;
    defaults to `BROptions()`.

# Returns

- `action_dist::Vector{<:Integer}` : Updated `action_dist`.
"""
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

play!(brd::SamplingBRD, action::Integer, action_dist::Vector{<:Integer},
      options::BROptions=BROptions()) =
    play!(Random.GLOBAL_RNG, brd, action, action_dist, options)


# play

"""
    play(rng, brd, init_action_dist, options; num_reps)

Return the action distribution after `num_reps` times iteration

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : AbstractBRD instance.
- `init_action_dist::Vector{<:Integer}` : The initial distribution of actions
    for each players.
- `options::BROptions` : Options for `best_response` method;
    defaults to `BROptions()`.
- `num_reps::Integer` : The number of iterations; defaults to 1.

# Returns

- `::Vector{<:Integer}` : The action distribution after iterations
"""
function play(rng::AbstractRNG,
              brd::AbstractBRD,
              init_action_dist::Vector{<:Integer},
              options::BROptions=BROptions();
              num_reps::Integer=1)
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

"""
    time_series(rng, brd, ts_length, init_action_dist, options)

Return the time series of action distribution.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `brd::AbstractBRD` : Instance of the model.
- `ts_length::Integer` : The length of time series.
- `init_action_dist::Vector{<:Integer}` : Initial action distribution.
- `options::BROptions` : Options for `best_response` method;
    defaults to `BROptions()`.

# Returns

- `::Matrix{<:Integer}` : Time series of action distribution.
"""
function time_series(rng::AbstractRNG,
                     brd::AbstractBRD,
                     ts_length::Integer,
                     init_action_dist::Vector{<:Integer},
                     options::BROptions=BROptions())
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

"""
    time_series(rng, brd, ts_length, init_actions, options)

Return the time series of action distribution. Initial action distribution is
choosed randomly.

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
    player_ind_seq = rand(rng, 1:brd.N, ts_length)
    nums_actions = ntuple(i -> brd.num_actions, brd.N)
    init_actions = random_pure_actions(nums_actions)
    action_dist = zeros(Int, brd.num_actions)
    for i in 1:brd.N
        action_dist[actions[i]] += 1
    end
    time_series(rng, brd, ts_length, action_dist, options)
end

time_series(brd::AbstractBRD, ts_length::Integer,
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, brd, ts_length, options)
