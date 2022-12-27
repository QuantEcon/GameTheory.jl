#=
Tools for normal form games.
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

    function Player{N,T}(payoff_array::Array{T,N}) where {N,T<:Real}
        if prod(size(payoff_array)) == 0
            throw(ArgumentError("every player must have at least one action"))
        end
        return new(payoff_array)
    end
end

Player(payoff_array::Array{T,N}) where {T<:Real,N} = Player{N,T}(payoff_array)

Player{N,T}(player::Player{N,S}) where {N,T,S} =
    Player(Array{T}(player.payoff_array))

Base.convert(::Type{T}, player::Player) where {T<:Player} =
    player isa T ? player : T(player)

"""
    Player(T, player)

Convert `player` into a new `Player` instance with eltype `T`.

# Arguments

- `T::Type`
- `player::Player`

# Returns

- `::Player` : `Player` instance with eltype `T`.
"""
Player(::Type{T}, player::Player{N}) where {T<:Real,N} = Player{N,T}(player)

Player(player::Player{N,T}) where {N,T} = Player{N,T}(player)

num_actions(p::Player) = size(p.payoff_array, 1)
num_opponents(::Player{N}) where {N} = N - 1

Base.summary(player::Player) =
    string(Base.dims2string(size(player.payoff_array)),
           " ",
           split(string(typeof(player)), ".")[end])

function Base.show(io::IO, player::Player)
    print(io, summary(player))
    println(io, ":")
    X = player.payoff_array
    if !haskey(io, :compact) && length(axes(X, 2)) > 1
        io = IOContext(io, :compact => true)
    end
    Base.print_array(io, X)
end

# delete_action

"""
    delete_action(player, action[, player_idx=1])

Return a new Player instance with the action(s) specified by `action` deleted
from the action set of the player specified by `player_idx`.

# Arguments

- `player::Player` : Player instance.
- `action::Union{PureAction,AbstractVector{<:PureAction}}`: The action(s) to be
  deleted.
- `player_idx::Integer` : Index of the player to delete action(s) for.

# Returns

- `::Player` : `Player` instance with the action(s) deleted as specified.

# Examples

```julia
julia> player = Player([3 0; 0 3; 1 1])
3×2 Player{2, Int64}:
 3  0
 0  3
 1  1

julia> delete_action(player, 3)
2×2 Player{2, Int64}:
 3  0
 0  3

julia> delete_action(player, 1, 2)
3×1 Player{2, Int64}:
 0
 3
 1
```
"""
function delete_action(player::Player{N,T}, action::AbstractVector{<:PureAction},
                       player_idx::Integer=1) where {N,T}
    sel = Any[Colon() for i in 1:N]
    sel[player_idx] = setdiff(axes(player.payoff_array, player_idx), action)
    payoff_array_new = player.payoff_array[sel...]::Array{T,N}
    return Player(payoff_array_new)
end

delete_action(player::Player, action::PureAction, player_idx::Integer=1) =
    delete_action(player, [action], player_idx)

# Reference
# https://stackoverflow.com/questions/53637511/how-to-delete-the-specific-row-of-n-dimensional-array-in-julia

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
- `opponents_actions::PureActionProfile` : Tuple of N-1 opponents' pure actions.

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

- `player::Player{2}` : Player instance.
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

- `player::Player{2}` : Player instance.
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

- `player::Player{1}` : Player instance.
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

Return true if `own_action` is a best response to `opponents_actions`.

# Arguments

- `player::Player` : Player instance.
- `own_action::PureAction` : Own pure action (integer).
- $(opponents_actions_docstring)
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
  false otherwise.
"""
function is_best_response(player::Player,
                          own_action::PureAction,
                          opponents_actions::Union{Action,ActionProfile,Nothing};
                          tol::Real=1e-8)
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
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool` : True if `own_action` is a best response to `opponents_actions`;
  false otherwise.
"""
function is_best_response(player::Player,
                          own_action::MixedAction,
                          opponents_actions::Union{Action,ActionProfile,Nothing};
                          tol::Real=1e-8)
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
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns

- `best_responses::Vector{Int}` : Vector containing all the best response
  actions.
"""
function best_responses(player::Player,
                        opponents_actions::Union{Action,ActionProfile,Nothing};
                        tol::Real=1e-8)
    payoffs = payoff_vector(player, opponents_actions)
    payoff_max = maximum(payoffs)
    best_responses = findall(x -> x >= payoff_max - tol, payoffs)
    return best_responses
end

"""
    best_response([rng=GLOBAL_RNG], player, opponents_actions;
                  tie_breaking=:smallest, tol=1e-8)

Return a best response action to `opponents_actions`.

# Arguments

- `rng::AbstractRNG=GLOBAL_RNG` : Random number generator; relevant only with
  `tie_breaking=:random`.
- `player::Player` : Player instance.
- $(opponents_actions_docstring)
- `tie_breaking::Symbol` : Control how to break a tie (see Returns for details).
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns

- `::Int` : If `tie_breaking=:smallest`, returns the best response action with
  the smallest index; if `tie_breaking=:random`, returns an action randomly
  chosen from the best response actions.
"""
function best_response(rng::AbstractRNG, player::Player,
                       opponents_actions::Union{Action,ActionProfile,Nothing};
                       tie_breaking::Symbol=:smallest,
                       tol::Real=1e-8)
    if tie_breaking == :smallest
        payoffs = payoff_vector(player, opponents_actions)
        return argmax(payoffs)
    elseif tie_breaking == :random
        brs = best_responses(player, opponents_actions; tol=tol)
        return rand(rng, brs)
    else
        throw(ArgumentError(
            "tie_breaking must be one of `:smallest` or `:random`"
        ))
    end
end

best_response(player::Player,
              opponents_actions::Union{Action,ActionProfile,Nothing};
              tie_breaking::Symbol=:smallest, tol::Real=1e-8) =
    best_response(Random.GLOBAL_RNG, player, opponents_actions,
                  tie_breaking=tie_breaking, tol=tol)

"""
    BROptions

Struct to contain options for `best_response`.

# Fields
- `tol::Real=1e-8` : Tolerance level.
- `tie_breaking::Symbol=:smallest` : `:smallest` or `:random`.
- `rng::AbstractRNG=GLOBAL_RNG` : Random number generator.
"""
@with_kw struct BROptions{T<:Real,TR<:AbstractRNG}
    tol::T = 1e-8
    tie_breaking::Symbol = :smallest
    rng::TR = Random.GLOBAL_RNG
end

"""
    best_response(player, opponents_actions, options)

Return a best response action to `opponents_actions` with options as specified
by a `BROptions` instance `options`.
"""
best_response(player::Player,
              opponents_actions::Union{Action,ActionProfile,Nothing},
              options::BROptions) =
    best_response(options.rng, player, opponents_actions;
                  tie_breaking=options.tie_breaking, tol=options.tol)

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
    players = Vector{Player{N,T}}(undef, N)
    for i in 1:N
        sz = ntuple(j -> nums_actions[i-1+j <= N ? i-1+j : i-1+j-N], N)
        players[i] = Player(zeros(T, sz))
    end
    return NormalFormGame(players)
end

NormalFormGame(nums_actions::NTuple{N,Int}) where {N} =
    NormalFormGame(Float64, nums_actions)

# Check that the shapes of the payoff arrays are consistent
function is_consistent(players::NTuple{N,Player{N,T}}) where {N,T}
    shape_1 = size(players[1].payoff_array)
    for i in 2:N
        shape = size(players[i].payoff_array)
        for j in 1:(N-i+1)
            shape[j] == shape_1[i-1+j] || return false
        end
        for j in (N-i+2):N
            shape[j] == shape_1[j-(N-i+1)] || return false
        end
    end
    return true
end

"""
    NormalFormGame(players)

Constructor of an N-player NormalFormGame with a tuple of N Player instances.

# Arguments

- `players::NTuple{N,Player}` : Tuple of Player instances.
"""
function NormalFormGame(players::NTuple{N,Player{N,T}}) where {N,T}
    is_consistent(players) ||
        throw(ArgumentError("shapes of payoff arrays must be consistent"))
    nums_actions = ntuple(i -> num_actions(players[i]), N)
    return NormalFormGame{N,T}(players, nums_actions)
end

"""
    NormalFormGame(players)

Constructor of an N-player NormalFormGame with a vector of N Player instances.

# Arguments

- `players::Vector{Player}` : Vector of Player instances.
"""
NormalFormGame(players::Vector{Player{N,T}}) where {N,T} =
    NormalFormGame(ntuple(i -> players[i], N))

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
    colons = Base.front(ntuple(j -> Colon(), M))

    size(payoffs)[end] != N && throw(ArgumentError(
        "length of the array in the last axis must be equal to
         the number of players"
    ))

    players = ntuple(
        i -> Player(permutedims(view(payoffs, colons..., i),
                                     (i:N..., 1:i-1...)::typeof(dims))
                   ),
        Val(N)
    )
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

"""
    NormalFormGame(payoffs)

Construct an N-player NormalFormGame with an N-dimensional array `payoffs` of
vectors, where `payoffs[a_1, a_2, ..., a_N]` contains a vector of N payoff
values, one for each player, for the action profile (a\\_1, a\\_2, ..., a\\_N).

# Arguments

- `payoffs::Array{AbstractVector{T<:Real}}` : Array with ndims=N containing
  payoff profiles as vectors.
"""
function NormalFormGame(payoffs::Array{TV,N}) where
    {T<:Real,N,TV<:AbstractVector{T}}
    g = NormalFormGame(T, size(payoffs))
    for (i, player) in enumerate(g.players)
        payoffs_permted = PermutedDimsArray(payoffs, (i:N..., 1:i-1...))
        for a in eachindex(player.payoff_array)
            player.payoff_array[a] = payoffs_permted[a][i]
        end
    end
    return g
end

function NormalFormGame{N,T}(g::NormalFormGame{N,S}) where {N,T,S}
    players_new = ntuple(i -> Player{N,T}(g.players[i]), Val(N))
    return NormalFormGame(players_new)
end

Base.convert(::Type{T}, g::NormalFormGame) where {T<:NormalFormGame} =
    g isa T ? g : T(g)

"""
    NormalFormGame(T, g)

Convert `g` into a new `NormalFormGame` instance with eltype `T`.

# Arguments

- `T::Type`
- `g::NormalFormGame`

# Returns

- `::NormalFormGame` : `NormalFormGame` instance with eltype `T`.
"""
NormalFormGame(::Type{T}, g::NormalFormGame{N}) where {T<:Real,N} =
    NormalFormGame{N,T}(g)

NormalFormGame(g::NormalFormGame{N,T}) where {N,T} = NormalFormGame{N,T}(g)

Base.summary(g::NormalFormGame) =
    string(Base.dims2string(g.nums_actions),
           " ",
           split(string(typeof(g)), ".")[end])

# payoff_profile_array

"""
    payoff_profile_array(g)

Return an N-dimensional array of vectors, whose (a\\_1, ..., a\\_N)-entry
contains a vector of N payoff values, one for each player, for the action profile
(a\\_1, ..., a\\_N).

# Arguments

- `g::NormalFormGame` : N-player `NormalFormGame` instance.

# Returns

- `::Array{Vector,N}` : Array of payoff profiles.
"""
function payoff_profile_array(g::NormalFormGame{N,T}) where {N,T}
    payoff_profile_array =
        map(index -> Vector{T}(undef, N), CartesianIndices(g.nums_actions))
    for (i, player) in enumerate(g.players)
        for index in CartesianIndices(g.nums_actions)
            payoff_profile_array[index][i] =
                player.payoff_array[(index.I[i:end]..., index.I[1:i-1]...)...]
        end
    end
    return payoff_profile_array
end

function Base.show(io::IO, g::NormalFormGame)
    print(io, summary(g))
end

function Base.print(io::IO, g::NormalFormGame)
    print(io, summary(g))
    println(io, ":")
    X = payoff_profile_array(g)
    if !haskey(io, :compact) && length(axes(X, 2)) > 1
        io = IOContext(io, :compact => true)
    end
    Base.print_array(io, X)
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
function Base.getindex(g::NormalFormGame{1}, index::Integer)
    return g.players[1].payoff_array[index]
end

function Base.setindex!(g::NormalFormGame{N},
                        payoff_profile::AbstractVector{<:Real},
                        index::Integer...) where N
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
function Base.setindex!(g::NormalFormGame{1},
                        payoff::Real,
                        index::Integer)
    g.players[1].payoff_array[index] = payoff
    return payoff
end

Base.setindex!(g::NormalFormGame{N},
               payoff_profile::NTuple{N}, index...) where N =
    setindex!(g, collect(payoff_profile), index...)

# Indexing with CartesianIndices
Base.getindex(g::NormalFormGame{N}, ci::CartesianIndex{N}) where {N} =
    g[to_indices(g, (ci,))...]
Base.setindex!(g::NormalFormGame{N}, v, ci::CartesianIndex{N}) where {N} =
    g[to_indices(g, (ci,))...] = v

# delete_action

"""
    delete_action(g, action, player_idx)

Return a new `NormalFormGame` instance with the action(s) specified by `action`
deleted from the action set of the player specified by `player_idx`.

# Arguments

- `g::NormalFormGame` : `NormalFormGame` instance.
- `action::Union{PureAction, AbstractVector{<:PureAction}}` : The action(s) to
  be deleted.
- `player_idx::Integer` : Index of the player to delete action(s) for.

# Returns

- `::NormalFormGame` : `NormalFormGame` instance with the action(s) deleted as
  specified.
"""
function delete_action(g::NormalFormGame{N},
                       action::AbstractVector{<:PureAction},
                       player_idx::Integer) where N
    players_new  = [delete_action(player, action,
                    player_idx-i+1>0 ? player_idx-i+1 : player_idx-i+1+N)
                    for (i, player) in enumerate(g.players)]
    return NormalFormGame(players_new)
end

delete_action(g::NormalFormGame, action::PureAction, player_idx::Integer) =
    delete_action(g, [action], player_idx)

# is_nash

function is_nash(g::NormalFormGame, action_profile::ActionProfile;
                 tol::Real=1e-8)
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
                 tol::Real=1e-8)
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
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool`
""" is_nash

# Trivial game with 1 player
"""
    is_nash(g, action; tol=1e-8)

Return true if `action` is a Nash equilibrium of a trivial game with 1 player.

# Arguments

- `g::NormalFormGame{1}` : Instance of 1-player NormalFormGame.
- `action::Action` : Integer (pure action) or vector of reals (mixed action).
- `tol::Float64` : Tolerance to be used to determine best response actions.

# Returns

- `::Bool`
"""
is_nash(g::NormalFormGame{1}, action::Action; tol::Real=1e-8) =
    is_best_response(g.players[1], action, nothing, tol=tol)

is_nash(g::NormalFormGame{1}, action_profile::ActionProfile;
        tol::Real=1e-8) = is_nash(g, action_profile..., tol=tol)

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

- `mixed_action::Vector{Float64}` : The mixed action representation of the given
  pure action.
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

- `g::NormalFormGame` : Instance of N-player NormalFormGame.
- `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

- `::Bool`
""" is_pareto_efficient

@doc """
    is_pareto_dominant(g, action_profile)

Return true if `action_profile` is Pareto dominant for game `g`.

# Arguments

- `g::NormalFormGame` : Instance of N-player NormalFormGame.
- `action_profile::PureActionProfile` : Tuple of N integers (pure actions).

# Returns

- `::Bool`
""" is_pareto_dominant

# is_dominated

"""
    is_dominated(player, action; tol=1e-8,
                 lp_solver=GameTheory.highs_optimizer_silent)

Determine whether `action` is strictly dominated by some mixed action.

# Arguments

- `player::Player` : Player instance.
- `action::PureAction` : Integer representing a pure action.
- `tol::Real` : Tolerance level used in determining domination.
- `lp_solver` : Linear programming solver to be used internally. Pass a
  `MathOptInterface.AbstractOptimizer` type (such as `HiGHS.Optimizer`) if no
  option is needed, or a function (such as `GameTheory.highs_optimizer_silent`)
  to supply options.

# Returns

- `::Bool` : True if `action` is strictly dominated by some mixed action; false
  otherwise.

"""
function is_dominated(
    ::Type{T}, player::Player, action::PureAction; tol::Real=1e-8,
    lp_solver=highs_optimizer_silent
) where {T<:Real}
    payoff_array = player.payoff_array
    m, n = size(payoff_array, 1) - 1, prod(size(payoff_array)[2:end])

    ind = trues(num_actions(player))
    ind[action] = false

    A_ub = Matrix{T}(undef, (m+1, n))  # transposed
    A_ub[1:end-1, :] .= reshape(selectdim(payoff_array, 1, action), (1, n))
    A_ub[1:end-1, :] -= reshape(selectdim(payoff_array, 1, ind), (m, n))
    A_ub[end, :] .= 1

    a_eq = ones(T, m+1)
    a_eq[end] = 0

    c = zeros(T, m+1)
    c[end] = 1

    CACHE = MOIU.UniversalFallback(MOIU.Model{T}())
    optimizer = MOIU.CachingOptimizer(CACHE, lp_solver())
    x = MOI.add_variables(optimizer, m+1)
    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
            MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}.(c, x), 0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    for j in 1:n
        MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction{T}(
                MOI.ScalarAffineTerm{T}.(A_ub[:, j], x), 0
            ),
            MOI.LessThan{T}(0)
        )
    end
    MOI.add_constraint(
        optimizer,
        MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}.(a_eq, x), 0),
        MOI.EqualTo{T}(1)
    )
    # Nonnegativity
    for i in 1:m
        a = zeros(T, m+1)
        a[i] = -1
        MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}.(a, x), 0),
            MOI.LessThan{T}(0)
        )
    end
    MOI.optimize!(optimizer)
    status = MOI.get(optimizer, MOI.TerminationStatus())

    if status == MOI.OPTIMAL
        return (MOI.get(optimizer, MOI.ObjectiveValue()) > tol)::Bool
    elseif status == MOI.INFEASIBLE
        return false
    else
        throw(ErrorException("Error: solution status $(status)"))
    end
end

function is_dominated(
    ::Type{T}, player::Player{1}, action::PureAction; tol::Real=1e-8,
    lp_solver=highs_optimizer_silent
) where {T<:Real}
        payoff_array = player.payoff_array
        return maximum(payoff_array) > payoff_array[action] + tol
end

is_dominated(
    player::Player, action::PureAction; tol::Real=1e-8,
    lp_solver=highs_optimizer_silent
) = is_dominated(Float64, player, action, tol=tol, lp_solver=lp_solver)

# dominated_actions

"""
    dominated_actions(player; tol=1e-8,
                      lp_solver=GameTheory.highs_optimizer_silent)

Return a vector of actions that are strictly dominated by some mixed actions.

# Arguments

- `player::Player` : Player instance.
- `tol::Real` : Tolerance level used in determining domination.
- `lp_solver::Union{Type{<:MathOptInterface.AbstractOptimizer},Function}` :
  Linear programming solver to be used internally. Pass a
  `MathOptInterface.AbstractOptimizer` type (such as `HiGHS.Optimizer`) if no
  option is needed, or a function (such as `GameTheory.highs_optimizer_silent`)
  to supply options.

# Returns

- `out::Vector{Int}` : Vector of integers representing pure actions, each
  of which is strictly dominated by some mixed action.

"""
function dominated_actions(
    ::Type{T}, player::Player; tol::Real=1e-8, lp_solver=highs_optimizer_silent
) where {T<:Real}
    out = Int[]
    for action = 1:num_actions(player)
        if is_dominated(T, player, action, tol=tol, lp_solver=lp_solver)
            append!(out, action);
        end
    end

    return out
end

dominated_actions(
    player::Player; tol::Real=1e-8, lp_solver=highs_optimizer_silent
) = dominated_actions(Float64, player, tol=tol, lp_solver=lp_solver)
