#=
Tools for fictitious play

=#


# AbstractGain #

"""
    AbstractGain

Abstract type representing the gain of fictitious play moodel.
"""
abstract type AbstractGain end

"""
    DecreasingGain

Type representing the gain decresing over time. Subtype of `AbstractGain`.
"""
struct DecreasingGain <: AbstractGain end

"""
    ConstantGain

Type representing the gain constant over time. Subtype of `AbstractGain`.
"""
struct ConstantGain{T<:Real} <: AbstractGain
    size::T
end

step_size(T::Type, gain::DecreasingGain, t::Integer) = one(T)/(t+1)
step_size(T::Type, gain::ConstantGain, t::Integer) = T(gain.size)


# AbstractFictitiousPlay #

"""
    AbstractFictitiousPlay

Abstract type representing the fictitious play model.
"""
abstract type AbstractFictitiousPlay{N,T<:Real} end

"""
    FictitiousPlay{N, T, TG}

Type representing the fictitious play model with N players.
Subtype of `AbstractFictitiousPlay`.

# Fields

- `players::NTuple{N,Player{N,T}}` : Tuple of player instances.
- `nums_actions::NTuple{N,Int}` : Tuple of integers which are the number of
    actions for each player.
- `gain::TG` : Type of gain.
"""
struct FictitiousPlay{N,T<:Real,TG<:AbstractGain} <: AbstractFictitiousPlay{N,T}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
    gain::TG
end

"""
    FictitiousPlay(g[, gain])

Construct a `FictitiousPlay` instance.

# Arguments

- `g::NormalFormGame` : `NormalFormGame` instance.
- `gain::AbstractGain` : Argument to specify the gain or step size;
    `DecreasingGain()` (default) or `ConstantGain(size)`.

# Returns

- `::FictitiousPlay`
"""
FictitiousPlay(g::NormalFormGame, gain::AbstractGain=DecreasingGain()) =
    FictitiousPlay(g.players, g.nums_actions, gain)

"""
    FictitiousPlay(fp[, gain])

Construct a new `FictitiousPlay` instance from `fp`.

# Arguments

- `g::AbstractFictitiousPlay`
- `gain::AbstractGain`

# Returns

- `::FictitiousPlay`
"""
FictitiousPlay(fp::AbstractFictitiousPlay, gain::AbstractGain=fp.gain) =
    FictitiousPlay(fp.players, fp.nums_actions, gain)

"""
    StochasticFictitiousPlay{N, T, TG, TD}

Type representing the stochastic fictitious play model with N players.
Subtype of `AbstractFictitiousPlay`.

# Fields

- `players::NTuple{N,Player{N,T}}` : Tuple of player instances.
- `nums_actions::NTuple{N,Int}` : Tuple of integers which are the number of
    actions for each player.
- `gain::TG` : Type of gain.
- `d::TD` Distribution of the payoff shocks.
"""
struct StochasticFictitiousPlay{N,T<:Real,TD<:Distribution,
                                TG<:AbstractGain} <: AbstractFictitiousPlay{N,T}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
    d::TD
    gain::TG
end

"""
    StochasticFictitiousPlay(g, d[, gain])

Construct a `StochasticFictitiousPlay` instance.

# Arguments

- `g::NormalFormGame` : `NormalFormGame` instance.
- `d::Distributions.Distribution` : `Distribution` instance from which payoff
    perturbations are drawn.
- `gain::AbstractGain` : Argument to specify the gain or step size;
    `DecreasingGain()` (default) or `ConstantGain(size)`.

# Returns

- `::StochasticFictitiousPlay`
"""
StochasticFictitiousPlay(g::NormalFormGame, d::Distribution,
                         gain::AbstractGain=DecreasingGain()) =
    StochasticFictitiousPlay(g.players, g.nums_actions, d, gain)

"""
    StochasticFictitiousPlay(fp[, d, gain])

Construct a new `StochasticFictitiousPlay` instance from `fp`.

# Arguments

- `g::AbstractFictitiousPlay`
- `d::Distributions.Distribution`
- `gain::AbstractGain`

# Returns

- `::StochasticFictitiousPlay`
"""
StochasticFictitiousPlay(fp::AbstractFictitiousPlay, d::Distribution,
                         gain::AbstractGain=fp.gain) =
    StochasticFictitiousPlay(fp.players, fp.nums_actions, d, gain)
StochasticFictitiousPlay(fp::StochasticFictitiousPlay,
                         gain::AbstractGain=fp.gain) =
    StochasticFictitiousPlay(fp.players, fp.nums_actions, fp.d, gain)


# play!

"""
    play!(rng, fp, actions, options, brs, t)

Update `actions` which represents the normalized action history for each player.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::FictitiousPlay{N}` : FictitiousPlay instance.
- `actions::MixedActionProfile{TA,N}` : Normalized action history for each
    player.
- `options::BROptions` : Options for `best_response` method.
- `brs::Vector{Int}` : Vector used temporarily.
- `t::Integer` : Integer representing period.

# Returns

- `actions::MixedActionProfile` : Updated `actions`.
"""
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

"""
    play!(rng, fp, actions, options, brs, t)

Update `actions` which represents the normalized action history for each player.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::StochasticFictitiousPlay{N}` : StochasticFictitiousPlay instance.
- `actions::MixedActionProfile{TA,N}` : Normalized action history for each player.
- `options::BROptions` : Options for `best_response` method.
- `brs::Vector{Int}` : Vector used temporarily.
- `t::Integer` : Integer representing period.

# Returns

- `actions::MixedActionProfile` : Updated `actions`.
"""
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

"""
    play!(rng, fp, actions, options, num_reps, t_init)

Update normalized action history `num_reps` times.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : StochasticFictitiousPlay instance.
- `actions::MixedActionProfile{TA,N}` : Normalized action history for each player.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iteration.
- `t_init::Integer` : The period iteration starts.

# Returns

- `actions::MixedActionProfile` : Updated `actions`.
"""
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


# play

"""
    play(rng, fp, actions, options, num_reps, t_init)

Return normalized action history after `num_reps` times iterations.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : StochasticFictitiousPlay instance.
- `actions::MixedActionProfile{TA,N}` : Normalized action history for each player.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iteration.
- `t_init::Integer` : The period iteration starts.

# Returns

- `actions::MixedActionProfile` : Normalized action history after iterations.
"""
function play(rng::AbstractRNG=Random.GLOBAL_RNG,
              fp::AbstractFictitiousPlay{N},
              actions::MixedActionProfile{TA,N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N,TA<:Real}
    Tout = typeof(zero(TA)/one(TA))
    actions_copied::NTuple{N,Vector{Tout}} =
        ntuple(i -> copyto!(similar(actions[i], Tout), actions[i]), N)
    play!(rng, fp, actions_copied, options, num_reps=num_reps, t_init=t_init)
end

"""
    play(rng, fp, actions, options, num_reps, t_init)

Return normalized action history after `num_reps` times iterations.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : StochasticFictitiousPlay instance.
- `actions::PureActionProfile{TA,N}` : Normalized action history for each player.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iteration.
- `t_init::Integer` : The period iteration starts.

# Returns

- `::MixedActionProfile` : Normalized action history after iterations.
"""
function play(rng::AbstractRNG=Random.GLOBAL_RNG,
              fp::AbstractFictitiousPlay{N},
              actions::PureActionProfile{N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N}
    mixed_actions = ntuple(i -> pure2mixed(fp.nums_actions[i], actions[i]), N)
    play!(rng, fp, mixed_actions, options, num_reps=num_reps, t_init=t_init)
end

"""
    play(rng, fp, actions, options, num_reps, t_init)

Return normalized action history after `num_reps` times iterations.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : StochasticFictitiousPlay instance.
- `options::BROptions` : Options for `best_response` method.
- `num_reps::Integer` : The number of iteration.
- `t_init::Integer` : The period iteration starts.

# Returns

- `::MixedActionProfile` : Normalized action history after iterations.
"""
function play(rng::AbstractRNG=Random.GLOBAL_RNG,
              fp::AbstractFictitiousPlay{N},
              options::BROptions=BROptions();
              num_reps::Integer=1, t_init::Integer=1) where {N}
    play!(rng, fp, random_mixed_actions(rng, fp.nums_actions), options,
          num_reps=num_reps, t_init=t_init)
end


# time_series!

"""
    time_series!(rng, fp, out, options, t_init)

Update the `out` which

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : AbstractFictitiousPlay instance.
- `out::NTuple{N,Matrix{<:Real}}` : Tuple of matrices which represent the time
    series of normalized action history.
- `options::BROptions` : Options for `best_response`.
- `t_init::Integer` : The period iteration starts.

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
    @eval function time_series(rng::AbstractRNG=Random.GLOBAL_RNG,
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

@doc """
    time_series(rng, fp, ts_length, init_actions, options, t_init)

Return time series of normalized action histories.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N,T}` : AbstractFictitiousPlay instance.
- `ts_length::Integer` : The length of time series.
- `init_actions::MixedActionProfile{TA,N}` : Initial action history.
- `options::BROptions` : Options for `best_response` method.
- `t_init::Integer` : The period iteration starts.

# Returns

- `::NTuple{N,Matrix{<:Real}}` : Time series of normalized action history.
"""

@doc """
    time_series(rng, fp, ts_length, init_actions, options, t_init)

Return time series of normalized action histories.

# Arguments

- `rng::AbstractRNG` : Random number generator.
- `fp::AbstractFictitiousPlay{N,T}` : AbstractFictitiousPlay instance.
- `ts_length::Integer` : The length of time series.
- `init_actions::PureActionProfile{TA,N}` : Initial action history.
- `options::BROptions` : Options for `best_response` method.
- `t_init::Integer` : The period iteration starts.

# Returns

- `::NTuple{N,Matrix{<:Real}}` : Time series of normalized action history.
"""

"""
    time_series(rng, fp, ts_length, options, t_init)

Return time series of normalized action histories.

# Arguments

- `rng::AbstractRNG` : Random number generator used.
- `fp::AbstractFictitiousPlay{N}` : AbstractFictitiousPlay instance.
- `ts_length::Integer` : The length of time series.
- `options::BROptions` : Options for `best_response` method.
- `t_init::Integer` : The period iteration starts.

# Returns
- `::NTuple{N,Matrix{<:Real}}` : Time series of normalized action history.
"""
function time_series(rng::AbstractRNG=Random.GLOBAL_RNG,
                     fp::AbstractFictitiousPlay{N},
                     ts_length::Integer,
                     options::BROptions=BROptions();
                     t_init::Integer=1) where {N}
    time_series(rng, fp, ts_length, random_mixed_actions(rng, fp.nums_actions),
                options, t_init=t_init)
end
