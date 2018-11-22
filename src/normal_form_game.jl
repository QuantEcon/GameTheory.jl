#=
Tools for normal form games.

Authors: Daisuke Oyama

=#

const opponents_actions_docstring = """
`opponents_actions::Union{Action,ActionProfile,Nothing}` : Profile of N-1
  opponents' actions. If N=2, then it must be a vector of reals (in which case
  it is treated as the opponent's mixed action) or a scalar of integer (in which
  case it is treated as the opponent's pure action). If N>2, then it must be a
  tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions).
  (For the degenerate case N=1, it must be `nothing`.)"""


# Player #

"""
    Player{N,T}

Type representing a player in an N-player normal form game.

# Fields

- `payoff_array::Array{T<:Real}` : Array representing the player's payoff
  function.
"""
struct Player{N,T<:Real}
    payoff_array::Array{T,N}
end

num_actions(p::Player) = size(p.payoff_array, 1)
num_opponents(::Player{N}) where {N} = N - 1

Base.summary(player::Player) =
    string(Base.dims2string(size(player.payoff_array)),
           " ",
           split(string(typeof(player)), ".")[end])

function Base.show(io::IO, player::Player)
    print(io, summary(player))
    println(io, ":")
    Base.print_array(io, player.payoff_array)
end

# payoff_vector

# To resolve definition ambiguity
function payoff_vector(player::Player, opponents_actions::Tuple{})
    throw(ArgumentError("input tuple must not be empty"))
end

"""
    payoff_vector(player, opponents_actions)

Return a vector of payoff values for a Player in an N>2 player game, one for
each own action, given a tuple of the opponents' pure actions.

# Arguments

- `player::Player` : Player instance.
- `opponents_actions::PureActionProfile` : Tuple of N-1 opponents' pure
  actions.

# Returns

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
    payoff_vector(player, opponents_actions)

Return a vector of payoff values for a Player in an N>2 player game, one for
each own action, given a tuple of the opponents' mixed actions.

# Arguments

- `player::Player` : Player instance.
- `opponents_actions::MixedActionProfile` : Tuple of N-1 opponents' mixed
  actions.

# Returns

- `::Vector` : Payoff vector.
"""
function payoff_vector(player::Player{N,T1},
                       opponents_actions::MixedActionProfile{T2}) where {N,T1,T2}
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
    payoff_vector(player, opponent_action)

Return a vector of payoff values for a Player in a 2-player game, one for each
own action, given the opponent's pure action.

# Arguments

- `player::Player` : Player instance.
- `opponent_action::PureAction` : Opponent's pure action (integer).

# Returns

- `::Vector` : Payoff vector.
"""
function payoff_vector(player::Player{2}, opponent_action::PureAction)
    return player.payoff_array[:, opponent_action]
end

"""
    payoff_vector(player, opponent_action)

Return a vector of payoff values for a Player in a 2-player game, one for each
own action, given the opponent's mixed action.

# Arguments

- `player::Player` : Player instance.
- `opponent_action::MixedAction` : Opponent's mixed action (vector of reals).

# Returns

- `::Vector` : Payoff vector.
"""
function payoff_vector(player::Player{2}, opponent_action::MixedAction)
    # player.num_opponents == 1
    return player.payoff_array * opponent_action
end

# Trivial case with player.num_opponents == 0
"""
    payoff_vector(player, opponent_action)

Return a vector of payoff values for a Player in a trivial game with 1 player,
one for each own action.

# Arguments

- `player::Player` : Player instance.
- `opponent_action::Nothing`

# Returns

- `::Vector` : Payoff vector.
"""
function payoff_vector(player::Player{1}, opponent_action::Nothing)
    return player.payoff_array
end

# _reduce_ith_opponent

# Given an N-d array `payoff_array` of shape (n_{1}, ..., n_{i+1}, 1, ..., 1),
# return the N-d array of payoffs of shape (n_{1}, ..., n_{i}, 1, 1, ..., 1)
# given by fixing the i-th opponent's pure or mixed action to be `action`.
for (S, ex_mat_action) in ((PureAction, :(A[:, action])),
                           (MixedAction, :(A * action)))
    @eval function _reduce_ith_opponent(payoff_array::Array{T,N},
                                        i::Int, action::$S) where {N,T}
        shape = size(payoff_array)
        A = reshape(payoff_array, (prod(shape[1:i]), shape[i+1]))::Matrix{T}
        out = $(ex_mat_action)
        shape_new = tuple(shape[1:i]..., ones(Int, N-i)...)::NTuple{N,Int}
        return reshape(out, shape_new)
    end
end

# is_best_response

"""
    is_best_response(player, own_action, opponents_actions; tol=1e-8)

Return True if `own_action` is a best response to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- `own_action::PureAction` : Own pure action (integer).
- $(opponents_actions_docstring)
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
  false otherwise.
"""
function is_best_response(player::Player,
                          own_action::PureAction,
                          opponents_actions::Union{Action,ActionProfile,Nothing};
                          tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    return payoffs[own_action] >= payoff_max - tol
end

"""
    is_best_response(player, own_action, opponents_actions; tol=1e-8)

Return true if `own_action` is a best response to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- `own_action::MixedAction` : Own mixed action (vector of reals).
- $(opponents_actions_docstring)
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
  false otherwise.
"""
function is_best_response(player::Player,
                          own_action::MixedAction,
                          opponents_actions::Union{Action,ActionProfile,Nothing};
                          tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    return dot(own_action, payoffs) >= payoff_max - tol
end

# best_response

"""
    best_responses(player, opponents_actions; tol=1e-8)

Return all the best response actions to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `best_responses::Vector{Int}` : Vector containing all the best response
  actions.
"""
function best_responses(player::Player,
                        opponents_actions::Union{Action,ActionProfile,Nothing};
                        tol::Float64=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    best_responses = findall(x -> x >= payoff_max - tol, payoffs)
    return best_responses
end

"""
    best_response(player, opponents_actions; tie_breaking="smallest", tol=1e-8)

Return a best response action to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `tie_breaking::AbstractString("smallest")` : Control how to break a tie (see
  Returns for details).
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Int` : If tie_breaking="smallest", returns the best response action with
  the smallest index; if tie_breaking="random", returns an action randomly
  chosen from the best response actions.
"""
function best_response(player::Player,
                       opponents_actions::Union{Action,ActionProfile,Nothing};
                       tie_breaking::AbstractString="smallest",
                       tol::Float64=1e-8)
    if tie_breaking == "smallest"
        payoffs = payoff_vector(player, opponents_actions)
        return argmax(payoffs)
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
    best_response(player, opponents_actions, payoff_perturbation)

Return the perturbed best response to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `payoff_perturbation::Vector{Float64}` : Vector of length equal to the number
  of actions of the player containing the values ("noises") to be added to the
  payoffs in determining the best response.

# Returns

- `::Int` : Best response action.
"""
function best_response(player::Player,
                       opponents_actions::Union{Action,ActionProfile,Nothing},
                       payoff_perturbation::Vector{Float64})
    length(payoff_perturbation) != num_actions(player) &&
        throw(ArgumentError(
            "length of payoff_perturbation must be $(num_actions(player))"
        ))

    payoffs = payoff_vector(player, opponents_actions) + payoff_perturbation
    return argmax(payoffs)
end


# NormalFormGame #

"""
    NormalFormGame{N,T}

Type representing an N-player normal form game.

# Fields

- `players::NTuple{N,Player{N,T<:Real}}` : Tuple of Player instances.
- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
"""
struct NormalFormGame{N,T<:Real}
    players::NTuple{N,Player{N,T}}
    nums_actions::NTuple{N,Int}
end

num_players(::NormalFormGame{N}) where {N} = N

function NormalFormGame(::Tuple{})  # To resolve definition ambiguity
    throw(ArgumentError("input tuple must not be empty"))
end

"""
    NormalFormGame([T=Float64], nums_actions)

Constructor of an N-player NormalFormGame, consisting of payoffs all 0.

# Arguments

- `T::Type` : Type of payoff values; defaults to `Float64` if not specified.
- `nums_actions::NTuple{N,Int}` : Numbers of actions of the N players.
"""
function NormalFormGame(T::Type, nums_actions::NTuple{N,Int}) where N
    # TODO: can we still get inference to work but avoid the `::NTuple` below?
    players::NTuple{N,Player{N,T}} =
        ntuple(i -> Player(zeros(tuple(nums_actions[i:end]...,
                                       nums_actions[1:i-1]...))),
               N)
    return NormalFormGame{N,T}(players, nums_actions)
end

NormalFormGame(nums_actions::NTuple{N,Int}) where {N} =
    NormalFormGame(Float64, nums_actions)

"""
    NormalFormGame(players)

Constructor of an N-player NormalFormGame with a tuple of N Player instances.

# Arguments

- `players::NTuple{N,Player}` : Tuple of Player instances.
"""
function NormalFormGame(players::NTuple{N,Player{N,T}}) where {N,T}
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
    NormalFormGame(players)

Constructor of an N-player NormalFormGame with a vector of N Player instances.

# Arguments

- `players::Vector{Player}` : Vector of Player instances.
"""
NormalFormGame(players::Vector{Player{N,T}}) where {N,T} =
    NormalFormGame(tuple(players...)::NTuple{N,Player{N,T}})

"""
    NormalFormGame(players...)

Constructor of an N-player NormalFormGame with N Player instances.

# Arguments

- `players::Player{N,T}...` : N Player instances

# Examples

```julia
# p1, p2, and p3 are all of type `Player{3,T}` for some `T`
NormalFormGame(p1, p2, p3)
```
"""
function NormalFormGame(players::Player{N,T}...) where {N,T}
    length(players) != N && error("Need $N players")
    NormalFormGame(players)  # use constructor for Tuple of players above
end

"""
    NormalFormGame(payoffs)

Construct an N-player NormalFormGame for N>=2 with an array `payoffs` of M=N+1
dimensions, where `payoffs[a_1, a_2, ..., a_N, :]` contains a profile of N
payoff values.

# Arguments

- `payoffs::Array{T<:Real}` : Array with ndims=N+1 containing payoff profiles.
"""
function NormalFormGame(payoffs::Array{T,M}) where {T<:Real,M}
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
    NormalFormGame(payoffs)

Construct a symmetric 2-player NormalFormGame with a square matrix.

# Arguments

- `payoffs::Matrix{T<:Real}` : Square matrix representing each player's payoff
  matrix.
"""
function NormalFormGame(payoffs::Matrix{T}) where T<:Real
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

function Base.getindex(g::NormalFormGame{N,T},
                       index::Integer...) where {N,T}
    length(index) != N &&
        throw(DimensionMismatch("index must be of length $N"))

    payoff_profile = Array{T}(undef, N)
    for (i, player) in enumerate(g.players)
        payoff_profile[i] =
            player.payoff_array[(index[i:end]..., index[1:i-1]...)...]
    end
    return payoff_profile
end

# Trivial game with 1 player
function Base.getindex(g::NormalFormGame{1,T}, index::Integer) where T
    return g.players[1].payoff_array[index]
end

function Base.setindex!(g::NormalFormGame{N,T},
                        payoff_profile::Vector{S},
                        index::Integer...) where {N,T,S<:Real}
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
function Base.setindex!(g::NormalFormGame{1,T},
                        payoff::S,
                        index::Integer) where {T,S<:Real}
    g.players[1].payoff_array[index] = payoff
    return payoff
end

# Indexing with CartesianIndices
Base.getindex(g::NormalFormGame{N}, ci::CartesianIndex{N}) where {N} =
    g[to_indices(g, (ci,))...]
Base.setindex!(g::NormalFormGame{N}, v, ci::CartesianIndex{N}) where {N} =
    g[to_indices(g, (ci,))...] = v

# is_nash

function is_nash(g::NormalFormGame, action_profile::ActionProfile;
                 tol::Float64=1e-8)
    for (i, player) in enumerate(g.players)
        own_action = action_profile[i]
        opponents_actions =
            tuple(action_profile[i+1:end]..., action_profile[1:i-1]...)
        if !(is_best_response(player, own_action, opponents_actions, tol=tol))
            return false
        end
    end
    return true
end

function is_nash(g::NormalFormGame{2}, action_profile::ActionProfile;
                 tol::Float64=1e-8)
    for (i, player) in enumerate(g.players)
        own_action, opponent_action =
            action_profile[i], action_profile[3-i]
        if !(is_best_response(player, own_action, opponent_action, tol=tol))
            return false
        end
    end
    return true
end

# attach docstring for N>1 player games
@doc """
    is_nash(g, action_profile; tol=1e-8)

Return true if `action_profile` is a Nash equilibrium.

# Arguments

- `g::NormalFormGame` : Instance of N-player NormalFormGame.
- `action_profile::ActionProfile` : Tuple of N integers (pure actions) or N
  vectors of reals (mixed actions).
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool`
""" is_nash

# Trivial game with 1 player
"""
    is_nash(g, action; tol=1e-8)

Return true if `action` is a Nash equilibrium of a trivial game with 1 player.

# Arguments

- `g::NormalFormGame` : Instance of 1-player NormalFormGame.
- `action::Action` : Integer (pure action) or vector of reals (mixed action).
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool`
"""
is_nash(g::NormalFormGame{1}, action::Action; tol::Float64=1e-8) =
    is_best_response(g.players[1], action, nothing, tol=tol)

is_nash(g::NormalFormGame{1}, action_profile::ActionProfile;
        tol::Float64=1e-8) = is_nash(g, action_profile..., tol=tol)

# Utility functions

"""
    pure2mixed(num_actions, action)

Convert a pure action to the corresponding mixed action.

# Arguments

- `num_actions::Integer` : The number of the pure actions (= the length of a
  mixed action).
- `action::PureAction` : The pure action to convert to the corresponding mixed
  action.

# Returns

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
        for profile in CartesianIndices(g.nums_actions)
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
    is_pareto_efficient(g, action_profile)

Return true if `action_profile` is Pareto efficient for game `g`.

# Arguments

* `g::NormalFormGame` : Instance of N-player NormalFormGame.
* `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

* `::Bool`
""" is_pareto_efficient

@doc """
    is_pareto_dominant(g, action_profile)

Return true if `action_profile` is Pareto dominant for game `g`.

# Arguments

* `g::NormalFormGame` : Instance of N-player NormalFormGame.
* `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

* `::Bool`
""" is_pareto_dominant

# is_dominated

"""
    is_dominated(player, action[, tol=1e-8, lp_solver=ClpSolver()])

Determine whether `action` is strictly dominated by some mixed
action.

# Arguments

- `player::Player` : Player instance.
- `action::PureAction` : Integer representing a pure action.
- `tol::Real` : Tolerance to be used.
- `lp_solver::AbstractMathProgSolver` : Allows users to choose a particular
  solver for linear programming problems. Options include ClpSolver(),
  CbcSolver(), GLPKSolverLP() and GurobiSolver(). By default, it choooses
  ClpSolver().

# Returns

- `::Bool` : True if `action` is strictly dominated by some mixed action;
  False otherwise.

"""
function is_dominated(player::Player{N,T}, action::PureAction;
                      tol::Real=1e-8,
                      lp_solver::MathProgBase.AbstractMathProgSolver=
                      ClpSolver()) where {N,T<:Real}
    payoff_array = player.payoff_array
    S = typeof(zero(T)/one(T))

    m, n = size(payoff_array, 1) - 1, prod(size(player.payoff_array)[2:end])

    ind = trues(num_actions(player))
    ind[action] = false

    A = Array{S}(undef, n+1, m+1)
    A[1:n, 1:m] = transpose(reshape(-selectdim(payoff_array, 1, ind) .+
                                    selectdim(payoff_array, 1, action:action),
                                    (m, n)))
    A[1:n, m+1] .= 1
    A[n+1, 1:m] .= 1
    A[n+1, m+1] = 0

    b = zeros(S, n+1)
    b[end] = 1

    c = zeros(S, m+1)
    c[end] = -1

    sense = Array{Char}(undef, n+1)
    sense[1:n] .= '<'
    sense[n+1] = '='

    res = linprog(c, A, sense, b, lp_solver)

    if res.status == :Optimal
        return res.sol[end] > tol
    elseif res.status == :Infeasible
        return false
    else
        throw(ErrorException("Error: solution status $(res.status)"))
    end
end

function is_dominated(player::Player{1}, action::PureAction;
                      tol::Real=1e-8,
                      lp_solver::MathProgBase.AbstractMathProgSolver=
                      ClpSolver())
        payoff_array = player.payoff_array
        return maximum(payoff_array) > payoff_array[action] + tol
end

# dominated_actions

"""
    dominated_actions(player[, tol=1e-8, lp_solver=ClpSolver()])

Return a vector of actions that are strictly dominated by some mixed actions.

# Arguments

- `player::Player` : Player instance.
- `tol::Real` : Tolerance level used in determining domination.

# Returns

- `out::Vector{Int}` : Vector of integers representing pure actions, each
  of which is strictly dominated by some mixed action.

"""
function dominated_actions(player::Player; tol::Real=1e-8,
                           lp_solver::MathProgBase.AbstractMathProgSolver=
                           ClpSolver())
    out = Vector{Int}(undef, 0)
    for action = 1:num_actions(player)
        if is_dominated(player, action, tol=tol, lp_solver=lp_solver)
            append!(out, action);
        end
    end

    return out
end
