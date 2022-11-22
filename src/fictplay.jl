#=
Tools for fictitious play

=#


# AbstractGain #

"""
    AbstractGain

Abstract type representing the gain in a fictitious play model.
"""
abstract type AbstractGain end

"""
    DecreasingGain

Type representing a decresing gain.
"""
struct DecreasingGain <: AbstractGain end

"""
    ConstantGain

Type representing a constant gain.
"""
struct ConstantGain{T<:Real} <: AbstractGain
    size::T
end

step_size(T::Type, gain::DecreasingGain, t::Integer) = one(T)/(t+1)
step_size(T::Type, gain::ConstantGain, t::Integer) = T(gain.size)


# AbstractFictitiousPlay #

"""
    AbstractFictitiousPlay

Abstract type representing a fictitious play model.
"""
abstract type AbstractFictitiousPlay{N,T<:Real} end

"""
    FictitiousPlay{N, T, TG}

Type representing a fictitious play model with N players.

# Fields

- `players::NTuple{N,Player{N,T}}` : Tuple of `Player` instances.
- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for
  each player.
- `gain::TG<:AbstractGain` : Gain type.
"""
struct FictitiousPlay{N,T<:Real,TG<:AbstractGain} <: AbstractFictitiousPlay{N,T}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
    gain::TG
end

"""
    FictitiousPlay(g[, gain=DecreasingGain()])

Construct a `FictitiousPlay` instance from `NormalFormGame`.

# Arguments

- `g::NormalFormGame` : `NormalFormGame` instance.
- `gain::AbstractGain` : Argument to specify the gain or step size;
  `DecreasingGain()` or `ConstantGain(size)`.

# Returns

- `::FictitiousPlay` : The fictitious play model.
"""
FictitiousPlay(g::NormalFormGame, gain::AbstractGain=DecreasingGain()) =
    FictitiousPlay(g.players, g.nums_actions, gain)

"""
    FictitiousPlay(fp[, gain=fp.gain])

Construct a new `FictitiousPlay` instance from `fp`.

# Arguments

- `fp::AbstractFictitiousPlay` : `AbstractFictitiousPlay` instance.
- `gain::AbstractGain` : Argument to specify the gain or step size.

# Returns

- `::FictitiousPlay` : The fictitious play model.
"""
FictitiousPlay(fp::AbstractFictitiousPlay, gain::AbstractGain=fp.gain) =
    FictitiousPlay(fp.players, fp.nums_actions, gain)

"""
    StochasticFictitiousPlay{N, T, TG, TD}

Type representing a stochastic fictitious play model with N players.

# Fields

- `players::NTuple{N,Player{N,T}}` : Tuple of `Player` instances.
- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for
  each player.
- `gain::TG<:AbstractGain` : Gain type.
- `d::TD<:Distributions.Distribution` : `Distribution` instance from which
  payoff perturbations are drawn.
"""
struct StochasticFictitiousPlay{N,T<:Real,TD<:Distribution,
                                TG<:AbstractGain} <: AbstractFictitiousPlay{N,T}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
    d::TD
    gain::TG
end

"""
    StochasticFictitiousPlay(g, d[, gain=DecreasingGain()])

Construct a `StochasticFictitiousPlay` instance.

# Arguments

- `g::NormalFormGame` : `NormalFormGame` instance.
- `d::Distributions.Distribution` : `Distribution` instance from which payoff
  perturbations are drawn.
- `gain::AbstractGain` : Argument to specify the gain or step size;
  `DecreasingGain()` or `ConstantGain(size)`.

# Returns

- `::StochasticFictitiousPlay` : The stochastic fictitious play model.
"""
StochasticFictitiousPlay(g::NormalFormGame, d::Distribution,
                         gain::AbstractGain=DecreasingGain()) =
    StochasticFictitiousPlay(g.players, g.nums_actions, d, gain)

"""
    StochasticFictitiousPlay(fp[, d=fp.d, gain=fp.gain])

Construct a new `StochasticFictitiousPlay` instance from `fp`.

# Arguments

- `fp::AbstractFictitiousPlay` : `AbstractFictitiousPlay` instance.
- `d::Distributions.Distribution` : `Distribution` instance from which payoff
  perturbations are drawn.
- `gain::AbstractGain` : Argument to specify the gain or step size.

# Returns

- `::StochasticFictitiousPlay` : The stochastic fictitious play model.
"""
StochasticFictitiousPlay(fp::AbstractFictitiousPlay, d::Distribution,
                         gain::AbstractGain=fp.gain) =
    StochasticFictitiousPlay(fp.players, fp.nums_actions, d, gain)
StochasticFictitiousPlay(fp::StochasticFictitiousPlay,
                         gain::AbstractGain=fp.gain) =
    StochasticFictitiousPlay(fp.players, fp.nums_actions, fp.d, gain)


# play!

function play!(rng::AbstractRNG,
               fp::FictitiousPlay{N},
               actions::MixedActionProfile{TA,N},
               options::BROptions,
               brs::Vector{Int}, t::Integer) where {N,TA<:Real}
    for i in 1:N
        opponents_actions =
            tuple(actions[i+1:end]..., actions[1:i-1]...)
        brs[i] = best_response(fp.players[i], opponents_actions, options)
    end

    for i in 1:N
        actions[i] .*= 1 - step_size(TA, fp.gain, t)
        actions[i][brs[i]] += step_size(TA, fp.gain, t)
    end

    return actions
end

function play!(rng::AbstractRNG,
               fp::StochasticFictitiousPlay{N},
               actions::MixedActionProfile{TA,N},
               options::BROptions,
               brs::Vector{Int}, t::Integer) where {N,TA<:Real}
    for i in 1:N
        opponents_actions =
            tuple(actions[i+1:end]..., actions[1:i-1]...)
        perturbations = rand(rng, fp.d, fp.nums_actions[i])
        brs[i] = best_response(fp.players[i], opponents_actions, perturbations)
    end

    for i in 1:N
        actions[i] .*= 1 - step_size(TA, fp.gain, t)
        actions[i][brs[i]] += step_size(TA, fp.gain, t)
    end

    return actions
end

play!(fp::StochasticFictitiousPlay{N}, actions::MixedActionProfile{TA,N},
      options::BROptions, brs::Vector{Int}, t::Integer) where {N,TA<:Real} =
    play!(Random.GLOBAL_RNG, fp, actions, options, brs, t)

function play!(rng::AbstractRNG,
               fp::AbstractFictitiousPlay{N},
               actions::MixedActionProfile{TA,N},
               options::BROptions=BROptions();
               num_reps::Integer=1, t_init::Integer=1) where {N,TA<:Real}
    brs = Vector{Int}(undef, N)
    for t in t_init:(t_init+num_reps-1)
        play!(rng, fp, actions, options, brs, t)
    end
    return actions
end

@doc """
    play!(rng, fp, actions[, options=BROptions(); num_reps=1, t_init=1])

Update action profile `num_reps` times.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : `AbstractFictitiousPlay` instance.
- `actions::MixedActionProfile{TA,N}` : Mixed action profile for each player.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iterations.
- `t_init::Integer` : The period when the iteration starts.

# Returns

- `actions::MixedActionProfile` : Updated mixed action profile.
"""


# play

function play(rng::AbstractRNG,
              fp::AbstractFictitiousPlay{N},
              actions::MixedActionProfile{TA,N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N,TA<:Real}
    Tout = typeof(zero(TA)/one(TA))
    actions_copied::NTuple{N,Vector{Tout}} =
        ntuple(i -> copyto!(similar(actions[i], Tout), actions[i]), N)
    play!(rng, fp, actions_copied, options, num_reps=num_reps, t_init=t_init)
end

play(fp::AbstractFictitiousPlay, actions::MixedActionProfile,
     options::BROptions=BROptions(); num_reps::Integer=1, t_init::Integer=1) =
    play(Random.GLOBAL_RNG, fp, actions, options,
         num_reps=num_reps, t_init=t_init)

function play(rng::AbstractRNG,
              fp::AbstractFictitiousPlay{N},
              actions::PureActionProfile{N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N}
    mixed_actions = ntuple(i -> pure2mixed(fp.nums_actions[i], actions[i]), N)
    play!(rng, fp, mixed_actions, options, num_reps=num_reps, t_init=t_init)
end

play(fp::AbstractFictitiousPlay, actions::PureActionProfile,
     options::BROptions=BROptions(); num_reps::Integer=1, t_init::Integer=1) =
    play(Random.GLOBAL_RNG, fp, actions, options,
         num_reps=num_reps, t_init=t_init)

function play(rng::AbstractRNG,
              fp::AbstractFictitiousPlay{N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N}
    play!(rng, fp, random_mixed_actions(rng, fp.nums_actions), options,
          num_reps=num_reps, t_init=t_init)
end

play(fp::AbstractFictitiousPlay, options::BROptions=BROptions();
     num_reps::Integer=1, t_init::Integer=1) =
    play(Random.GLOBAL_RNG, fp, options, num_reps=num_reps, t_init=t_init)

@doc """
    play([rng=Random.GLOBAL_RNG, ]fp, actions[, options=BROptions();
          num_reps=1, t_init=1])

Return a new action profile after `num_reps` times iterations.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `actions::ActionProfile` : Action profile used in the initial period;
  `PureActionProfile`, `MixedActionProfile`, or nothing. If nothing, mixed
  action profile is randomly selected.
- `fp::AbstractFictitiousPlay{N}` : `AbstractFictitiousPlay` instance.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of the iterations.
- `t_init::Integer` : The period when the iteration starts.

# Returns

- `::MixedActionProfile` : The new action profile after iterations.
"""

# time_series!

"""
    time_series!(rng, fp, out[, options=BROptions(); t_init=1])

Update the tuple of matrices `out` which is used in `time_series` method.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : `AbstractFictitiousPlay` instance.
- `out::NTuple{N,Matrix{<:Real}}` : Tuple of matrices which represent the time
    series of mixed action profile.
- `options::BROptions` : Options for `best_response`.
- `t_init::Integer` : The period when the iteration starts.

# Returns

- `out::NTuple{N,Matrix{<:Real}}` : Updated `out`.
"""
function time_series!(rng::AbstractRNG,
                      fp::AbstractFictitiousPlay{N},
                      out::NTuple{N,Matrix{<:Real}},
                      options::BROptions=BROptions();
                      t_init::Integer=1) where {N}
    ts_length = size(out[1], 2)
    actions = ntuple(i -> out[i][:, 1], N)
    brs = Vector{Int}(undef, N)

    for j in 2:ts_length
        play!(rng, fp, actions, options, brs, t_init - 1 + j - 1)
        for i in 1:N
            out[i][:, j] = actions[i]
        end
    end

    return out
end


# time_series

function _copy_action_to!(dest::AbstractVector, src::MixedAction)
    dest[:] = src
    return dest
end

function _copy_action_to!(dest::AbstractVector, src::PureAction)
    dest .= 0
    dest[src] = 1
    return dest
end

for (ex_TAS, ex_where, ex_T) in (
        (:(MixedActionProfile{TA,N}), (:(N), :(T<:Real), :(TA<:Real)), :(TA)),
        (:(PureActionProfile{N}), (:(N), :(T<:Real)), :(T))
    )
    @eval function time_series(rng::AbstractRNG,
                               fp::AbstractFictitiousPlay{N,T},
                               ts_length::Integer,
                               init_actions::$(ex_TAS),
                               options::BROptions=BROptions();
                               t_init::Integer=1) where $(ex_where...)
        Tout = typeof(zero($(ex_T))/one($(ex_T)))
        out::NTuple{N,Matrix{Tout}} =
            ntuple(i -> Matrix{Tout}(undef, fp.nums_actions[i], ts_length), N)
        for i in 1:N
            _copy_action_to!(@views(out[i][:, 1]), init_actions[i])
        end
        time_series!(rng, fp, out, options, t_init=t_init)
    end
end

time_series(fp::AbstractFictitiousPlay, ts_length::Integer,
            init_actions::ActionProfile, options::BROptions=BROptions();
            t_init::Integer=1) =
    time_series(Random.GLOBAL_RNG, fp, ts_length, init_actions, options,
                t_init=t_init)

function time_series(rng::AbstractRNG,
                     fp::AbstractFictitiousPlay{N},
                     ts_length::Integer,
                     options::BROptions=BROptions();
                     t_init::Integer=1) where {N}
    time_series(rng, fp, ts_length, random_mixed_actions(rng, fp.nums_actions),
                options, t_init=t_init)
end

time_series(fp::AbstractFictitiousPlay, ts_length::Integer,
            options::BROptions=BROptions(); t_init::Integer=1) =
    time_series(Random.GLOBAL_RNG, fp, ts_length, options, t_init=t_init)

@doc """
    time_series([rng=Random.GLOBAL_RNG, ]fp, ts_length, init_actions
                [, options=BROptions(); t_init=1])

Return a time series of mixed action profiles.

# Arguments

- `rng::AbstractRNG` : Random number generator.
- `fp::AbstractFictitiousPlay{N,T}` : `AbstractFictitiousPlay` instance.
- `ts_length::Integer` : The length of the time series.
- `init_actions::ActionProfile` : Action profile used in the initial period;
  `PureActionProfile`, `MixedActionProfile`, or nothing. If nothing, mixed
  action profile is randomly selected.
- `options::BROptions` : Options for `best_response` method.
- `t_init::Integer` : The period when the iteration starts.

# Returns

- `::NTuple{N,Matrix{<:Real}}` : The time series of mixed action profiles.
"""
