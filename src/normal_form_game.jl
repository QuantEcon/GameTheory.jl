#=
Tools for normal form games.

Authors: Daisuke Oyama

=#

const opponents_actions_docstring = """
`opponents_actions::Union{Action,ActionProfile,Void}` : Profile of N-1
opponents' actions. If N=2, then it must be a vector of reals (in which case
it is treated as the opponent's mixed action) or a scalar of integer (in which
case it is treated as the opponent's pure action). If N>2, then it must be a
tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions).
(For the degenerate case N=1, it must be `nothing`.)"""


# Player #

"""
Type representing a player in an N-player normal form game.

##### Arguments

- `payoff_array::Array{T<:Real}` : Array representing the player's payoff
function.

##### Fields

- `payoff_array::Array{T<:Real}` : Array representing the player's payoff
function.

"""
struct Player{N,T<:Real}
    payoff_array::Array{T,N}
end

num_actions(p::Player) = size(p.payoff_array, 1)
num_opponents{N}(::Player{N}) = N - 1

Base.summary(player::Player) =
    string(Base.dims2string(size(player.payoff_array)),
           " ",
           split(string(typeof(player)), ".")[end])

function Base.show(io::IO, player::Player)
    print(io, summary(player))
    println(io, ":")
    Base.showarray(io, player.payoff_array, false, header=false)
end

# payoff_vector

# To resolve definition ambiguity
function payoff_vector(player::Player, opponents_actions::Tuple{})
    throw(ArgumentError("input tuple must not be empty"))
end

"""
Return a vector of payoff values for a Player in an N>2 player game, one for
each own action, given a tuple of the opponents' pure actions.

##### Arguments

- `player::Player` : Player instance.
- `opponents_actions::PureActionProfile` : Tuple of N-1 opponents' pure
actions.

##### Returns

- `::Vector` : Payoff vector.

"""
function payoff_vector(player::Player, opponents_actions::PureActionProfile)
    length(opponents_actions) != num_opponents(player) &&
        throw(ArgumentError(
            "length of opponents_actions must be $(num_opponents(player))"
        ))
    payoffs = player.payoff_array
    for i in num_opponents(player):-1:1
        payoffs = _reduce_ith_opponent(payoffs, i, opponents_actions[i])
    end
    return vec(payoffs)
end

"""
Return a vector of payoff values for a Player in an N>2 player game, one for
each own action, given a tuple of the opponents' mixed actions.

##### Arguments

- `player::Player` : Player instance.
- `opponents_actions::MixedActionProfile` : Tuple of N-1 opponents' mixed
actions.

##### Returns

- `::Vector` : Payoff vector.

"""
function payoff_vector{N,T1,T2}(player::Player{N,T1},
                                opponents_actions::MixedActionProfile{T2})
    length(opponents_actions) != num_opponents(player) &&
        throw(ArgumentError(
            "length of opponents_actions must be $(num_opponents(player))"
        ))
    S = promote_type(T1, T2)
    payoffs::Array{S,N} = player.payoff_array
    for i in num_opponents(player):-1:1
        payoffs = _reduce_ith_opponent(payoffs, i, opponents_actions[i])
    end
    return vec(payoffs)
end

"""
Return a vector of payoff values for a Player in a 2-player game, one for each
own action, given the opponent's pure action.

##### Arguments

- `player::Player` : Player instance.
- `opponent_action::PureAction` : Opponent's pure action (integer).

##### Returns

- `::Vector` : Payoff vector.

"""
function payoff_vector(player::Player{2}, opponent_action::PureAction)
    return player.payoff_array[:, opponent_action]
end

"""
Return a vector of payoff values for a Player in a 2-player game, one for each
own action, given the opponent's mixed action.

##### Arguments

- `player::Player` : Player instance.
- `opponent_action::MixedAction` : Opponent's mixed action (vector of reals).

##### Returns

- `::Vector` : Payoff vector.

"""
function payoff_vector(player::Player{2}, opponent_action::MixedAction)
    # player.num_opponents == 1
    return player.payoff_array * opponent_action
end

# Trivial case with player.num_opponents == 0
"""
Return a vector of payoff values for a Player in a trivial game with 1 player,
one for each own action.

##### Arguments

- `player::Player` : Player instance.
- `opponent_action::Void`

##### Returns

- `::Vector` : Payoff vector.

"""
function payoff_vector(player::Player{1}, opponent_action::Void)
    return player.payoff_array
end

# _reduce_ith_opponent

# Given an N-d array `payoff_array` of shape (n_{1}, ..., n_{i+1}, 1, ..., 1),
# return the N-d array of payoffs of shape (n_{1}, ..., n_{i}, 1, 1, ..., 1)
# given by fixing the i-th opponent's pure or mixed action to be `action`.
for (S, ex_mat_action) in ((PureAction, :(A[:, action])),
                           (MixedAction, :(A * action)))
    @eval function _reduce_ith_opponent{N,T}(payoff_array::Array{T,N},
                                             i::Int, action::$S)
        shape = size(payoff_array)
        A = reshape(payoff_array, (prod(shape[1:i]), shape[i+1]))
        out = $(ex_mat_action)
        shape_new = tuple(shape[1:i]..., ones(Int, N-i)...)::NTuple{N,Int}
        return reshape(out, shape_new)
    end
end

# is_best_response

"""
Return True if `own_action` is a best response to `opponents_actions`.

##### Arguments

- `player::Player` : Player instance.
- `own_action::PureAction` : Own pure action (integer).
- $(opponents_actions_docstring)
- `;tol::Float64` : Tolerance to be used to determine best response actions.

##### Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
valse otherwise.

"""
function is_best_response(player::Player,
                          own_action::PureAction,
                          opponents_actions::Union{Action,ActionProfile,Void};
                          tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    return payoffs[own_action] >= payoff_max - tol
end

"""
Return true if `own_action` is a best response to `opponents_actions`.

##### Arguments

- `player::Player` : Player instance.
- `own_action::MixedAction` : Own mixed action (vector of reals).
- $(opponents_actions_docstring)
- `;tol::Float64` : Tolerance to be used to determine best response actions.

##### Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
false otherwise.

"""
function is_best_response(player::Player,
                          own_action::MixedAction,
                          opponents_actions::Union{Action,ActionProfile,Void};
                          tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    return dot(own_action, payoffs) >= payoff_max - tol
end

# best_response

"""
Return all the best response actions to `opponents_actions`.

##### Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `;tol::Float64` : Tolerance to be used to determine best response actions.

##### Returns

- `best_responses::Vector{Int}` : Vector containing all the best response
actions.

"""
function best_responses(player::Player,
                        opponents_actions::Union{Action,ActionProfile,Void};
                        tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    best_responses = find(x -> x >= payoff_max - tol, payoffs)
    return best_responses
end

"""
Return a best response action to `opponents_actions`.

##### Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `tie_breaking::AbstractString("smallest")` : Control how to break a tie (see
Returns for details).
- `tol::Float64` : Tolerance to be used to determine best response actions.

##### Returns

- `::Int` : If tie_breaking="smallest", returns the best response action with
the smallest index; if tie_breaking="random", returns an action randomly chosen
from the best response actions.

"""
function best_response(player::Player,
                       opponents_actions::Union{Action,ActionProfile,Void};
                       tie_breaking::AbstractString="smallest",
                       tol::Float64=1e-8)
    if tie_breaking == "smallest"
        payoffs = payoff_vector(player, opponents_actions)
        return indmax(payoffs)
    elseif tie_breaking == "random"
        brs = best_responses(player, opponents_actions; tol=tol)
        return rand(brs)
    else
        throw(ArgumentError(
            "tie_breaking must be one of 'smallest' or 'random'"
        ))
    end
end

# Perturbed best response
"""
Return the perturbed best response to `opponents_actions`.

##### Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `payoff_perturbation::Vector{Float64}` : Vector of length equal to the number
of actions of the player containing the values ("noises") to be added to the
payoffs in determining the best response.

##### Returns

- `::Int` : Best response action.

"""
function best_response(player::Player,
                       opponents_actions::Union{Action,ActionProfile,Void},
                       payoff_perturbation::Vector{Float64})
    length(payoff_perturbation) != num_actions(player) &&
        throw(ArgumentError(
            "length of payoff_perturbation must be $(num_actions(player))"
        ))

    payoffs = payoff_vector(player, opponents_actions) + payoff_perturbation
    return indmax(payoffs)
end


# NormalFormGame #

"""
Class representing an N-player normal form game.

##### Fields

- `players::NTuple{N,Player{N,T<:Real}}` : Tuple of Player instances.
- `N::Int` : The number of players.
- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
player.

"""
struct NormalFormGame{N,T<:Real}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
end

num_players{N}(::NormalFormGame{N}) = N

function NormalFormGame(::Tuple{})  # To resolve definition ambiguity
    throw(ArgumentError("input tuple must not be empty"))
end

"""
Constructor of an N-player NormalFormGame, consisting of payoffs all 0.

##### Arguments

- `T::Type` : Type of payoff values; defaults to `Float64` if not specified.
- `nums_actions::NTuple{N,Int}` : Numbers of actions of the N players.

"""
function NormalFormGame{N}(T::Type, nums_actions::NTuple{N,Int})
    # TODO: can we still get inference to work but avoid the `::NTuple` below?
    players::NTuple{N,Player{N,T}} =
        ntuple(i -> Player(zeros(tuple(nums_actions[i:end]...,
                                       nums_actions[1:i-1]...))),
               N)
    return NormalFormGame{N,T}(players, nums_actions)
end

NormalFormGame{N}(nums_actions::NTuple{N,Int}) =
    NormalFormGame(Float64, nums_actions)

"""
Constructor of an N-player NormalFormGame.

##### Arguments

- `players::NTuple{N,Player}` : Tuple of Player instances.

"""
function NormalFormGame{N,T}(players::NTuple{N,Player{N,T}})
    # Check that the shapes of the payoff arrays are consistent
    shape_1 = size(players[1].payoff_array)
    for i in 2:N
        shape = size(players[i].payoff_array)
        if shape != tuple(shape_1[i:end]..., shape_1[1:i-1]...)
            throw(ArgumentError("shapes of payoff arrays must be consistent"))
        end
    end

    nums_actions::NTuple{N,Int} =
        tuple([num_actions(player) for player in players]...)
    return NormalFormGame{N,T}(players, nums_actions)
end

"""
Constructor of an N-player NormalFormGame.

##### Arguments

- `players::Vector{Player}` : Vector of Player instances.

"""
NormalFormGame{N,T}(players::Vector{Player{N,T}}) =
    NormalFormGame(tuple(players...)::NTuple{N,Player{N,T}})

"""
Constructor of an N-player NormalFormGame.

##### Arguments

- `players::Player{N,T}...` : N Player instances

##### Examples

```julia
# p1, p2, and p3 are all of type `Player{3,T}` for some `T`
NormalFormGame(p1, p2, p3)
```
"""
function NormalFormGame{N,T}(players::Player{N,T}...)
    length(players) != N && error("Need $N players")
    NormalFormGame(players)  # use constructor for Tuple of players above
end

"""
    NormalFormGame{T<:Real,M}(payoffs::Array{T,M})

Construct an N-player NormalFormGame for N>=2 with an array `payoffs` of M=N+1
dimensions, where `payoffs[a_1, a_2, ..., a_N, :]` contains a profile of N
payoff values.

# Arguments

* `payoffs::Array{T<:Real}` : Array with ndims=N+1 containing payoff profiles.
"""
function NormalFormGame{T<:Real,M}(payoffs::Array{T,M})
    N = M - 1
    dims = Base.front(size(payoffs))
    colons = Base.front(ntuple(j -> Colon(), M)::NTuple{M,Colon})

    size(payoffs)[end] != N && throw(ArgumentError(
        "length of the array in the last axis must be equal to
         the number of players"
    ))

    players = [
        Player(permutedims(view(payoffs, colons..., i),
                           (i:N..., 1:i-1...)::typeof(dims))
        ) for i in 1:N
    ]
    # Call NormalFormGame{N,T}(players::Vector{Player{N,T}})
    NormalFormGame(players)
end

"""
    NormalFormGame{T<:Real}(payoffs::Matrix{T})

Construct a symmetric 2-player NormalFormGame with a square matrix.

# Arguments

* `payoffs::Matrix{T<:Real}` : Square matrix representing each player's payoff
  matrix.
"""
function NormalFormGame{T<:Real}(payoffs::Matrix{T})
    n, m = size(payoffs)
    n != m && throw(ArgumentError(
        "symmetric two-player game must be represented by a square matrix"
    ))
    player = Player(payoffs)
    return NormalFormGame(player, player)
end

Base.summary(g::NormalFormGame) =
    string(Base.dims2string(g.nums_actions),
           " ",
           split(string(typeof(g)), ".")[end])

# TODO: add printout of payoff arrays
function Base.show(io::IO, g::NormalFormGame)
    print(io, summary(g))
end

function Base.getindex{N,T}(g::NormalFormGame{N,T},
                            index::Integer...)
    length(index) != N &&
        throw(DimensionMismatch("index must be of length $N"))

    payoff_profile = Array{T}(N)
    for (i, player) in enumerate(g.players)
        payoff_profile[i] =
            player.payoff_array[(index[i:end]..., index[1:i-1]...)...]
    end
    return payoff_profile
end

# Trivial game with 1 player
function Base.getindex{T}(g::NormalFormGame{1,T}, index::Integer)
    return g.players[1].payoff_array[index]
end

function Base.setindex!{N,T,S<:Real}(g::NormalFormGame{N,T},
                                     payoff_profile::Vector{S},
                                     index::Integer...)
    length(index) != N &&
        throw(DimensionMismatch("index must be of length $N"))
    length(payoff_profile) != N &&
        throw(DimensionMismatch("assignment must be of $N-array"))

    for (i, player) in enumerate(g.players)
        player.payoff_array[(index[i:end]...,
                             index[1:i-1]...)...] = payoff_profile[i]
    end
    return payoff_profile
end

# Trivial game with 1 player
function Base.setindex!{T,S<:Real}(g::NormalFormGame{1,T},
                                   payoff::S,
                                   index::Integer)
    g.players[1].payoff_array[index] = payoff
    return payoff
end

# Indexing with CartesianIndices
Base.getindex{N}(g::NormalFormGame{N}, ci::CartesianIndex{N}) =
    g[to_indices(g, (ci,))...]
Base.setindex!{N}(g::NormalFormGame{N}, v, ci::CartesianIndex{N}) =
    g[to_indices(g, (ci,))...] = v

# is_nash

function is_nash(g::NormalFormGame, action_profile::ActionProfile)
    for (i, player) in enumerate(g.players)
        own_action = action_profile[i]
        opponents_actions =
            tuple(action_profile[i+1:end]..., action_profile[1:i-1]...)
        if !(is_best_response(player, own_action, opponents_actions))
            return false
        end
    end
    return true
end

function is_nash(g::NormalFormGame{2}, action_profile::ActionProfile)
    for (i, player) in enumerate(g.players)
        own_action, opponent_action =
            action_profile[i], action_profile[3-i]
        if !(is_best_response(player, own_action, opponent_action))
            return false
        end
    end
    return true
end

# attach docstring for N>1 player games
@doc """
Return true if `action_profile` is a Nash equilibrium.

##### Arguments

- `g::NormalFormGame` : Instance of N-player NormalFormGame.
- `action_profile::ActionProfile` : Tuple of N integers (pure actions) or N
vectors of reals (mixed actions).

##### Returns

- `::Bool`

""" is_nash

# Trivial game with 1 player
"""
Return true if `action` is a Nash equilibrium of a trivial game with 1 player.

##### Arguments

- `g::NormalFormGame` : Instance of 1-player NormalFormGame.
- `action::Action` : Integer (pure action) or vector of reals (mixed action).

##### Returns

- `::Bool`

"""
is_nash(g::NormalFormGame{1}, action::Action) =
    is_best_response(g.players[1], action, nothing)

# Utility functions

"""
Convert a pure action to the corresponding mixed action.

##### Arguments

- `num_actions::Integer` : The number of the pure actions (= the length of a
mixed action).
- `action::PureAction` : The pure action to convert to the corresponding mixed
action.

##### Returns

- `mixed_action::Vector{Float64}` : The mixed action representation of the
given pure action.

"""
function pure2mixed(num_actions::Integer, action::PureAction)
    mixed_action = zeros(num_actions)
    mixed_action[action] = 1
    return mixed_action
end

# is_pareto_efficient & is_pareto_dominant

function pareto_inferior_to(payoff_profile1, payoff_profile2)
    all(payoff_profile2 .>= payoff_profile1) &&
    any(payoff_profile2 .> payoff_profile1)
end

function not_pareto_superior_to(payoff_profile1, payoff_profile2)
    any(payoff_profile2 .> payoff_profile1) ||
    all(payoff_profile2 .== payoff_profile1)
end

for (f, op) = ((:is_pareto_efficient, pareto_inferior_to), 
               (:is_pareto_dominant, not_pareto_superior_to))
    @eval function $(f)(g::NormalFormGame,
                        action_profile::PureActionProfile)
        payoff_profile0 = g[action_profile...]
        for profile in CartesianRange(g.nums_actions)
            if CartesianIndex(action_profile) != profile
                if ($(op)(payoff_profile0, g[profile]))
                    return false
                end
            end
        end
        return true
    end
end

@doc """
    is_pareto_efficient(g::NormalFormGame, action_profile::PureActionProfile)

Return true if `action_profile` is Pareto efficient for game `g`.

# Arguments

* `g::NormalFormGame` : Instance of N-player NormalFormGame.
* `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

* `::Bool`
""" is_pareto_efficient

@doc """
    is_pareto_dominant(g::NormalFormGame, action_profile::PureActionProfile)

Return true if `action_profile` is Pareto dominant for game `g`.

# Arguments

* `g::NormalFormGame` : Instance of N-player NormalFormGame.
* `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

* `::Bool`
""" is_pareto_dominant
