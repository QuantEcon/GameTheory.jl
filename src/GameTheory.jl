module GameTheory

# stdlib
using LinearAlgebra, Random
using Markdown

# Packages
using QuantEcon
using Combinatorics
using Parameters
using Distributions

# Optimization packages
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
using Clp

# Geometry packages
using Polyhedra
using CDDLib
using LRSLib

# Type aliases #

"""
    PureAction

Alias for `Integer`.
"""
const PureAction = Integer

"""
    MixedAction{T}

Alias for `Vector{T}` where `T<:Real`.
"""
const MixedAction{T<:Real} = Vector{T}

"""
    Action{T}

Alias for `Union{PureAction,MixedAction{T}}` where `T<:Real`.
"""
const Action{T<:Real} = Union{PureAction,MixedAction{T}}

"""
    PureActionProfile{N,T}

Alias for `NTuple{N,T}` where `T<:PureAction`.
"""
const PureActionProfile{N,T<:PureAction} = NTuple{N,T}

"""
    MixedActionProfile{N,T}

Alias for `NTuple{N,MixedAction{T}}` where `T<:Real`.
"""
const MixedActionProfile{N,T<:Real} = NTuple{N,MixedAction{T}}

"""
    ActionProfile{N,T}

Alias for `Union{PureActionProfile{N,T},MixedActionProfile{N,T}}`.
"""
const ActionProfile{N,T} = Union{PureActionProfile{N,T},MixedActionProfile{N,T}}

const RatOrInt = Union{Rational,Int}

# package code goes here
include("normal_form_game.jl")
include("homotopy_continuation.jl")
include("lrsnash.jl")
include("pure_nash.jl")
include("repeated_game.jl")
include("random.jl")
include("lemke_howson.jl")
include("support_enumeration.jl")
include("vertex_enumeration.jl")
include("util.jl")
include("generators/Generators.jl")

include("fictplay.jl")
include("localint.jl")
include("brd.jl")
include("logitdyn.jl")

export
    # Types
    Player, NormalFormGame,

    # Type aliases
    Action, MixedAction, PureAction,
    ActionProfile, MixedActionProfile, PureActionProfile,

    # Normal form game functions
    best_response, best_responses, is_best_response, payoff_vector,
    is_nash, pure2mixed, pure_strategy_NE, is_pareto_efficient,
    is_pareto_dominant, is_dominated, dominated_actions, delete_action,
    payoff_profile_array,

    # General functions
    num_players, num_actions, num_opponents,

    # Utilities
    BROptions,

    # Nash Equilibrium
    pure_nash,

    # Repeated Games
    RepeatedGame, unpack, flow_u_1, flow_u_2, flow_u, best_dev_i,
    best_dev_1, best_dev_2, best_dev_payoff_i, best_dev_payoff_1,
    best_dev_payoff_2, worst_value_i, worst_value_1, worst_value_2,
    worst_values, outerapproximation, AS, uniquetolrows,

    # Random Games
    random_game, covariance_game,
    random_pure_actions, random_mixed_actions,

    # Lemke-Howson
    lemke_howson,

    # Support Enumeration
    support_enumeration, support_enumeration_task,

    # Vertex Enumeration
    vertex_enumeration, vertex_enumeration_task,

    # LRS
    lrsnash,

    # Homotopy Continuation
    hc_solve,

    # Learning algorithms
    play!, play, time_series,
    AbstractGain, DecreasingGain, ConstantGain,
    AbstractFictitiousPlay, FictitiousPlay, StochasticFictitiousPlay,
    AbstractRevision, SimultaneousRevision, AsynchronousRevision,
    LocalInteraction,
    AbstractBRD, BRD, KMR, SamplingBRD,
    LogitDynamics

end # module
