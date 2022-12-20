#= 
Tools for local interaction model

=#

using SparseArrays


# AbstractRevision #

"""
    AbstractRevision

Abstract type representing revision method.
"""
abstract type AbstractRevision end

"""
    AsynchronousRevision

Type representing an asynchronous revision.
"""
struct AsynchronousRevision <: AbstractRevision end

"""
    SimultaneousRevision

Type representing a simultaneous revision.
"""
struct SimultaneousRevision <: AbstractRevision end


# LocalInteraction

"""
    LocalInteraction{N, T, S, A, TR}

Type representing the local interaction model with N players.

# Fields

- `players::NTuple{N,Player{2,T}}` : Tuple of `Player` instances.
- `num_actions::Integer` : The number of actions for players.
- `adj_matrix::Array{S,2}` : Adjacency matrix of the graph in the model.
- `revision<:AbstractRevision` : The way to revise the action profile.
"""
struct LocalInteraction{N,T<:Real,S<:Real,A<:Integer,TR<:AbstractRevision}
    players::NTuple{N,Player{2,T}}
    num_actions::Int
    adj_matrix::SparseMatrixCSC{S,A}
    revision::TR
end

"""
    LocalInteraction(g, adj_matrix[, revision=SimultaneousRevision()])

Construct a `LocalInteraction` instance.

# Arguments

- `g::NormalFormGame` : The game used in the model.
- `adj_matrix::AbstractMatrix` : Adjacency matrix of the graph in the model.
- `revision::AbstractRevision` : Arguments to specify the revision method;
  `SimultaneousRevision()` or `AsynchronousRevision()`.

# Returns

- `::LocalInteraction` : The local interaction model.
"""
function LocalInteraction(g::NormalFormGame{2,T},
                          adj_matrix::AbstractMatrix{S},
                          revision::AbstractRevision=SimultaneousRevision()
                          ) where {T<:Real,S<:Real}
    if size(adj_matrix, 1) != size(adj_matrix, 2)
        throw(ArgumentError("Adjacency matrix must be square"))
    end
    N = size(adj_matrix, 1)
    players = ntuple(i -> g.players[1], N)
    num_actions = g.nums_actions[1]
    if num_actions != g.nums_actions[2]
        throw(ArgumentError("Payoff matrix must be square"))
    end
    sparse_adj = sparse(adj_matrix)::SparseMatrixCSC{S}
    return LocalInteraction(players, num_actions, sparse_adj, revision)
end

"""
    LocalInteraction(payoff_matrix,
                     adj_matrix[, revision=SimultaneousRevision()])

Construct a `LocalInteraction` instance.

# Arguments

- `payoff_matrix::Matrix` : The payoff matrix of the game.
- `adj_matrix::AbstractMatrix` : Adjacency matrix of the graph in the model.
- `revision::AbstractRevision` : Arguments to specify the revision method.
  `SimultaneousRevision()` or `AsynchronousRevision`

# Returns

- `::LocalInteraction` : The local interaction model.
"""
function LocalInteraction(payoff_matrix::Matrix{T},
                          adj_matrix::AbstractMatrix{S},
                          revision::AbstractRevision=SimultaneousRevision()
                          ) where {T<:Real,S<:Real}
    N = size(adj_matrix, 1)
    if N != size(adj_matrix, 2)
        throw(ArgumentError("Adjacency matrix must be square"))
    end
    players = ntuple(i -> Player(payoff_matrix), N)
    num_actions = size(payoff_matrix, 1)
    if num_actions != size(payoff_matrix, 2)
        throw(ArgumentError("Payoff matrix must be square"))
    end
    sparse_adj = sparse(adj_matrix)::SparseMatrixCSC{S}
    return LocalInteraction(players, num_actions, sparse_adj, revision)
end


# play!

function play!(li::LocalInteraction{N},
               actions::Vector{<:Integer},
               player_ind::AbstractVector{<:Integer},
               options::BROptions) where N
    actions_matrix = sparse(1:N, actions, ones(Int, N), N, li.num_actions)
    opponent_action = li.adj_matrix[player_ind,:] * actions_matrix
    for (k, i) in enumerate(player_ind)
        actions[i] = best_response(li.players[i], Vector(opponent_action[k,:]),
                                   options)
    end
    return actions
end

play!(li::LocalInteraction, actions::Vector{<:Integer}, player_ind::Integer,
      options::BROptions) = play!(li, actions, [player_ind], options)

play!(li::LocalInteraction{N}, actions::Vector{<:Integer},
      options::BROptions) where {N} = play!(li, actions, 1:N, options)

function play!(li::LocalInteraction,
               actions::Vector{<:Integer},
               options::BROptions,
               player_ind::Union{AbstractVector{<:Integer},Integer},
               num_reps::Integer=1)
    for t in 1:num_reps
        play!(li, actions, player_ind, options)
    end
    return actions
end

@doc """
    play!(li, actions, options, player_ind[, num_reps=1])

Update an action profile `num_reps` times.

# Arguments

- `li::LocalInteraction` : `LocalInteraction` instance.
- `actions::Vector{<:Integer}` : Action profile in the intial period.
- `options::BROptions` : Options for `best_response` method.
- `player_ind::Union{Vector{<:Integer},Integer} : Integer or vector of integers
  representing the index of players to take an action.
- `num_reps::Integer` : The number of iterations.

# Returns

- `actions::Vector{Int}` : Updated `actions`.
"""


# play

function play(li::LocalInteraction{N},
              actions::PureActionProfile,
              player_ind::Union{AbstractVector{<:Integer},Integer},
              options::BROptions=BROptions();
              num_reps::Integer=1) where N
    actions_vector = [i for i in actions]
    actions_vector = play!(li, actions_vector, options, player_ind, num_reps)
    new_actions = ntuple(i -> actions_vector[i], N)
    return new_actions
end

function play(li::LocalInteraction{N},
              actions::PureActionProfile,
              options::BROptions=BROptions();
              num_reps::Integer=1) where N
    play(li, actions, 1:N, options, num_reps=num_reps)
end

@doc """
    play(li, actions, player_ind[, options=BROptions(); num_reps=1])

Return the action profile after `num_reps` time iterations.

# Arguments

- `li::LocalInteraction{N}` : `LocalInteraction` instance.
- `actions::PureActionProfile` : Initial actions of each players.
- `player_ind::Union{AbstractVector{<:Integer},Integer}` : Integer or vector of
  integers representing the index of players to take an action with asynchronous
  revision.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iterations.

# Returns

- `::PureActionProfile` : Actions of each players after iterations.
"""


# time_series!

"""
    time_series!(li, out, options, player_ind_seq)

Update the matrix `out` which is used in `time_series` method given player index
sequence.

# Arguments

- `li::LocalInteraction{N}` : `LocalInteraction` instance.
- `out::Matrix{<:Integer}` : Matrix representing a time series of action
  profiles.
- `options::BROptions` : Options for `best_response` method.
- `player_ind_seq::Vector{<:Integer}` : Vector representing the index of players
  to take an action.

# Returns

- `out::Matrix{<:Integer}` : Updated `out`.
"""
function time_series!(li::LocalInteraction{N},
                      out::Matrix{<:Integer},
                      options::BROptions,
                      player_ind_seq::Vector{<:Integer}) where N
    ts_length = size(out, 2)
    if ts_length != length(player_ind_seq) + 1
        throw(ArgumentError("The length of `ts_length` and
                             `player_ind_seq` are mismatched"))
    end

    actions = [out[i,1] for i in 1:N]
    for t in 2:ts_length
        play!(li, actions, options, player_ind_seq[t-1])
        out[:,t] = actions
    end

    return out
end

"""
    time_series!(li, out, options)

Update the matrix `out` which is used in `time_series` method. All players take
their actions simultaneously.

# Arguments

- `li::LocalInteraction{N}` : `LocalInteraction` instance.
- `out::Matrix{<:Integer}` : Matrix representing a time series of action
  profiles.
- `options::BROptions` : Options for `best_response` method.

# Returns

- `out::Matrix{<:Integer}` : Updated `out`.
"""
function time_series!(li::LocalInteraction{N},
                      out::Matrix{<:Integer},
                      options::BROptions) where N
    ts_length = size(out, 2)
    actions = [out[i,1] for i in 1:N]
    for t in 2:ts_length
        play!(li, actions, options)
        out[:,t] = actions
    end
    
    return out
end


# time_series

function time_series(rng::AbstractRNG,
                     li::LocalInteraction{N},
                     ts_length::Integer,
                     init_actions::PureActionProfile,
                     player_ind_seq::Vector{<:Integer},
                     options::BROptions=BROptions()) where N
    out = Matrix{Int}(undef, N, ts_length)
    for i in 1:N
        out[i,1] = init_actions[i]
    end
    time_series!(li, out, options, player_ind_seq)
end

time_series(li::LocalInteraction, ts_length::Integer,
            init_actions::PureActionProfile, player_ind_seq::Vector{<:Integer},
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, li, ts_length, init_actions, player_ind_seq,
                options)

function time_series(rng::AbstractRNG,
                     li::LocalInteraction{N,T,S,A,TR},
                     ts_length::Integer,
                     init_actions::PureActionProfile,
                     options::BROptions=BROptions()
                    ) where {N,T,S,A,TR<:SimultaneousRevision}
    out = Matrix{Int}(undef, N, ts_length)
    for i in 1:N
        out[i, 1] = init_actions[i]
    end
    time_series!(li, out, options)
end

time_series(li::LocalInteraction{N,T,S,A,TR},
            ts_length::Integer,
            init_actions::PureActionProfile,
            options::BROptions=BROptions()
           ) where {N,T,S,A,TR<:SimultaneousRevision} =
    time_series(Random.GLOBAL_RNG, li, ts_length, init_actions, options)

function time_series(rng::AbstractRNG,
                     li::LocalInteraction{N,T,S,A,TR},
                     ts_length::Integer,
                     init_actions::PureActionProfile,
                     options::BROptions=BROptions()
                    ) where {N,T,S,A,TR<:AsynchronousRevision}
    player_ind_seq = rand(rng, 1:N, ts_length-1)
    time_series(rng, li, ts_length, init_actions, player_ind_seq, options)
end

time_series(li::LocalInteraction{N,T,S,A,TR},
            ts_length::Integer,
            init_actions::PureActionProfile,
            options::BROptions=BROptions()
           ) where {N,T,S,A,TR<:AsynchronousRevision} =
    time_series(Random.GLOBAL_RNG, li, ts_length, init_actions, options)

function time_series(rng::AbstractRNG,
                     li::LocalInteraction{N},
                     ts_length::Integer,
                     player_ind_seq::Vector{<:Integer},
                     options::BROptions=BROptions()) where N
    nums_actions = ntuple(i -> li.num_actions, N)
    init_actions = random_pure_actions(rng, nums_actions)
    time_series(rng, li, ts_length, init_actions, player_ind_seq, options)
end

time_series(li::LocalInteraction, ts_length::Integer,
            player_ind_seq::Vector{<:Integer}, options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, li, ts_length, player_ind_seq, options)

function time_series(rng::AbstractRNG,
                     li::LocalInteraction{N},
                     ts_length::Integer,
                     options::BROptions=BROptions()) where N
    nums_actions = ntuple(i -> li.num_actions, N)
    init_actions = random_pure_actions(rng, nums_actions)
    time_series(rng, li, ts_length, init_actions, options)
end

time_series(li::LocalInteraction, ts_length::Integer,
            options::BROptions=BROptions()) =
    time_series(Random.GLOBAL_RNG, li, ts_length, options)

@doc """
    time_series([rng=Random.GLOBAL_RNG, ]li, ts_length, init_actions,
                player_ind_seq[, options=BROptions()])

Return the time series of action profiles.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `li::LocalInteraction{N}` : `LocalInteraction` instance.
- `ts_length::Integer` : The length of time series.
- `init_actions::PureActionProfile` : Action profile in the initial period. If
  not provided, it is selected randomly.
- `player_ind_seq::Vector{<:Integer}` : Vector of integers representing the
  index of players to take an action with asynchronous revision. If not
  provided, it is selected randomly.
- `options::BROptions` : Options for `best_response` method.

# Returns

- `::Matrix{<:Integer}` : The time series of action profiles.
"""
