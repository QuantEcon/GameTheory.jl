module GameTheory

# stdlib
using LinearAlgebra, Random

# Packages
using QuantEcon
using Combinatorics
using Parameters

# Optimization packages
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
using HiGHS

# Geometry packages
using Polyhedra

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
MixedAction{T<:Real} = Vector{T}

"""
    Action{T}

Alias for `Union{PureAction,MixedAction{T}}` where `T<:Real`.
"""
Action{T<:Real} = Union{PureAction,MixedAction{T}}

"""
    PureActionProfile{N,T}

Alias for `NTuple{N,T}` where `T<:PureAction`.
"""
PureActionProfile{N,T<:PureAction} = NTuple{N,T}

"""
    MixedActionProfile{T,N}

Alias for `NTuple{N,MixedAction{T}}` where `T<:Real`.
"""
MixedActionProfile{T<:Real,N} = NTuple{N,MixedAction{T}}

"""
    ActionProfile

Alias for `Union{PureActionProfile,MixedActionProfile}`.
"""
const ActionProfile = Union{PureActionProfile,MixedActionProfile}

const RatOrInt = Union{Rational,Int}

# package code goes here
include("normal_form_game.jl")
include("lrsnash.jl")
include("pure_nash.jl")
include("repeated_game.jl")
include("random.jl")
include("support_enumeration.jl")
include("vertex_enumeration.jl")
include("util.jl")
include("generators/Generators.jl")

export
    # Types
    Player, NormalFormGame,

    # Type aliases
    Action, MixedAction, PureAction, ActionProfile,

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
    worst_values, outerapproximation,

    # Random Games
    random_game, covariance_game,
    random_pure_actions, random_mixed_actions,

    # Support Enumeration
    support_enumeration, support_enumeration_task,

    # LRS
    lrsnash,

    # Vertex Enumeration
    vertex_enumeration

end # module
