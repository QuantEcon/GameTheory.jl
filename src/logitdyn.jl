#=
    Tools for Logit Response Dynamics
=#

# LogitDynamics #

"""
    LogitDynamics{N, T, S}

Type representing the Logit-Dynamics model.

# Fields

- `players::NTuple{N,Player{N,T}}` : Tuple of `Player` instances.
- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `beta<:Real` : The level of noise in a player's decision.
- `choice_probs::Vector{Array}` : The choice probabilities of each action, one
  for each player.
"""
struct LogitDynamics{N,T<:Real,S<:Real}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
    beta::S
    choice_probs::Vector{Array}
end

"""
    LogitDynamics(g, beta)

Construct a `LogitDynamics` instance.

# Arguments

- `g::NormalFormGame{N,T}` : `NormalFormGame` instance.
- `beta::S` : The level of noise in players' decision.

# Returns

- `::LogitDynamics` : The Logit-Dynamics model.
"""
function LogitDynamics(g::NormalFormGame{N,T}, beta::S) where {N,T<:Real,S<:Real}
    choice_probs = Vector{Array}(undef, N)
    for (i, player) in enumerate(g.players)
        payoff_array = permutedims(player.payoff_array, vcat(2:N, 1))
        payoff_array_normalized = payoff_array .- maximum(payoff_array, dims=N)
        choice_probs[i] = cumsum(exp.(payoff_array_normalized .* beta),
                                 dims=N)
    end
    return LogitDynamics(g.players, g.nums_actions, beta, choice_probs)
end

"""
    play!(rng, ld, player_ind, actions)

Return a new action of player indexed by `player_ind` given each players' choice
probabilities.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `ld::LogitDynamics{N}` : `LogitDynamics` instance.
- `player_ind::Integer` : A player index who takes an action.
- `actions::Vector{<:Integer}` : The action profile.

# Returns

- `::Integer` : The new action of the player indexed by `player_ind`.
"""
function play!(rng::AbstractRNG, ld::LogitDynamics{N}, player_ind::Integer,
               actions::Vector{<:Integer}) where N
    oppponent_actions = [actions[player_ind+1:N]..., actions[1:player_ind-1]...]
    cdf = ld.choice_probs[player_ind][oppponent_actions..., :]
    random_value = rand(rng)
    next_action = searchsortedfirst(cdf, random_value*cdf[end])
    return next_action
end

"""
    play([rng=Random.GLOBAL_RNG,] ld, init_actions[; num_reps=1])

Return new action profile after `num_reps` iterations.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `ld::LogitDynamics{N}` : `LogitDynamics` instance.
- `init_actions::PureActionProfile` : Initial action profile.
- `num_reps::Integer` : The number of iterations.

# Returns

- `::Vector{<:Integer}` : New action profile.
"""
function play(rng::AbstractRNG,
              ld::LogitDynamics{N},
              init_actions::PureActionProfile;
              num_reps::Integer=1) where N
    actions = [m for m in init_actions]
    player_ind_seq = rand(rng, 1:N, num_reps)
    for player_ind in player_ind_seq
        actions[player_ind] = play!(rng, ld, player_ind, actions)
    end
    return actions
end

play(ld::LogitDynamics, init_actions::PureActionProfile;
     num_reps::Integer=1) =
    play(Random.GLOBAL_RNG, ld, init_actions, num_reps=num_reps)

"""
    time_series!(rng, ld, out, player_ind_seq)

Update the matrix `out` which is used in `time_series` method given a player
index sequence.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `ld::LogitDynamics{N}` : `LogitDynamics` instance.
- `out::Matrix{<:Integer}` : Matrix representing the time series of action
  profiles.
- `player_ind_seq::Vector{<:Integer}` : The sequence of player index, which is
  determined randomly.

# Returns

- `::Matrix{<:Integer}` : Updated `out`.
"""
function time_series!(rng::AbstractRNG,
                      ld::LogitDynamics{N},
                      out::Matrix{<:Integer},
                      player_ind_seq::Vector{<:Integer}) where N
    ts_length = size(out, 2)
    current_actions = [out[i, 1] for i in 1:N]
    for t in 1:ts_length-1
        current_actions[player_ind_seq[t]] = play!(rng, ld, player_ind_seq[t],
                                                   current_actions)
        for i in 1:N
            out[i, t+1] = current_actions[i]
        end
    end
    return out
end

"""
    time_series([rng=Random.GLOBAL_RNG,] ld, ts_length, init_actions)

Return a time series of action profiles.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `ld::LogitDynamics{N}` : `LogitDynamics` instance.
- `ts_length::Integer` : The length of time series.
- `init_actions::PureActionProfile` : Initial action profile.

# Returns

- `::Matrix{<:Integer}` : The time series of action profiles.
"""
function time_series(rng::AbstractRNG,
                     ld::LogitDynamics{N},
                     ts_length::Integer,
                     init_actions::PureActionProfile) where N
    player_ind_seq = rand(rng, 1:N, ts_length-1)
    out = Matrix{Int}(undef, N, ts_length)
    for i in 1:N
        out[i, 1] = init_actions[i]
    end
    time_series!(rng, ld, out, player_ind_seq)
end

time_series(ld::LogitDynamics, ts_length::Integer,
            init_actions::PureActionProfile) =
    time_series(Random.GLOBAL_RNG, ld, ts_length, init_actions)
