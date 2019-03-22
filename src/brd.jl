#=
    Tools for best response dynamics

=#

using StatsBase

abstract type AbstractBRD{T<:Real} end

struct BRD{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
end

BRD(N::Integer, payoff_array::Matrix{T}) where {T<:Real} =
    BRD(N, Player(payoff_array), size(payoff_array, 1))

struct KMR{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
    epsilon::Float64
end

KMR(N::Integer, payoff_array::Matrix{T},
    epsilon::Float64) where {T<:Real} =
    KMR(N, Player(payoff_array), size(payoff_array, 1), epsilon)

struct SamplingBRD{T<:Real} <: AbstractBRD{T}
    N::Int
    player::Player{2,T}
    num_actions::Int
    k::Int  #sample size
end

SamplingBRD(N::Integer, payoff_array::Matrix{T}, k::Integer) where {T<:Real} =
    SamplingBRD(N, Player(payoff_array), size(payoff_array, 1), k)

function set_action_dist(brd::AbstractBRD, actions::PureActionProfile)
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

function set_action_dist(brd::AbstractBRD)
    nums_actions = ntuple(i -> brd.num_actions, brd.N)
    actions = random_pure_actions(nums_actions)
    return set_action_dist(brd, actions)
end

function play!(brd::BRD, action::Integer, action_dist::Vector{<:Integer}, options::BROptions)
    action_dist[action] -= 1
    next_action = best_response(brd.player, actions, options)
    action_dist[next_action] += 1
    return action_dist
end

function play!(brd::KMR, action::Integer, action_dist::Vector{<:Integer}, options::BROptions)
    action_dist[action] -= 1
    if rand() <= brd.epsilon
        next_action = rand(1:brd.num_actions)
    else
        next_action = best_response(brd.player, actions, options)
    end
    action_dist[next_action] += 1
    return action_dist
end

function play!(brd::SamplingBRD, action::Integer, action_dist::Vector{<:Integer}, options::BROptions)
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

function time_series!(brd::AbstractBRD,
                      out::Matrix{<:Integer},
                      options::BROptions)
    ts_length = size(out, 1)
    player_ind_seq = rand(1:brd.N, ts_length)
    action_dist = [out[i,1] for i in 1:brd.num_actions]
    for t in 1:ts_length
        action = searchsortedlast(accumulate(+, action_dist), player_ind_seq[t])
        action_dist = play!(brd, action, action_dist, options)
        for i in 1:brd.num_actions
            out[i,t+1] = action_dist[i]
        end
    end
    return out
end

function time_series(brd::AbstractBRD,
                     ts_length::Integer,
                     init_actions::Union{PureActionProfile,Nothing}=nothing,
                     options::BROptions=BROptions())
    player_ind_seq = rand(1:brd.N, ts_length)
    action_dist = set_action_dist(brd, init_actions)
    out = Matrix{Int}(undef, brd.num_actions, ts_length)
    for i in 1:brd.num_actions
        out[i,1] = action_dist[i]
    end
    time_series!(brd, out, options)
end