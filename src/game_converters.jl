# PayoffVector #

"""
    PayoffLayout

Abstract supertype of the singleton types that specify the ordering of the
payoff values in a [`PayoffVector`](@ref).
"""
abstract type PayoffLayout end

"""
    PlayerMajor <: PayoffLayout

Player-major ordering, as in the GameTracer .gam format: all the payoffs of
player 1 come first, then those of player 2, and so on. Within each block,
action profiles are ordered with player 1 varying fastest, then player 2, ...,
player N (i.e., column-major order).
"""
struct PlayerMajor <: PayoffLayout end

"""
    ProfileMajor <: PayoffLayout

Profile-major ordering, as in the Gambit .nfg format: the payoffs of players
1, ..., N at the first action profile come first, then those at the second
action profile, and so on. Action profiles are ordered with player 1 varying
fastest, then player 2, ..., player N (i.e., column-major order).
"""
struct ProfileMajor <: PayoffLayout end

"""
    PayoffVector{L,N,T}

Internal intermediate representation that stores the payoffs of an `N`-player
game in a single flat vector of eltype `T`, ordered according to the layout
`L<:PayoffLayout`. See [`GAMPayoffVector`](@ref) and
[`NFGPayoffVector`](@ref) for the two layouts available.

Viewed as a `prod(nums_actions) × N` matrix whose `[a, i]` entry is the payoff
to player `i` at the `a`-th action profile in column-major order, `payoffs` is
that matrix vectorized in column-major order for `PlayerMajor`, and in
row-major order for `ProfileMajor`.

# Fields

- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `payoffs::Vector{T}` : Vector storing payoffs in the order specified by `L`.
"""
struct PayoffVector{L<:PayoffLayout,N,T<:Real}
    nums_actions::NTuple{N,Int}
    payoffs::Vector{T}

    function PayoffVector{L,N,T}(
        nums_actions::NTuple{N,Int}, payoffs::Vector{T}
    ) where {L<:PayoffLayout,N,T<:Real}
        any(n -> n <= 0, nums_actions) &&
            throw(ArgumentError("all nums_actions must be positive"))
        expected = prod(nums_actions) * N
        length(payoffs) == expected || throw(ArgumentError(
            "payoffs length mismatch: expected $expected, got $(length(payoffs))"
        ))
        return new(nums_actions, payoffs)
    end
end

# NOTE: The parameter bounds of an alias must be identical to those of
# `PayoffVector` (`T<:Real`), so that e.g. `GAMPayoffVector === PayoffVector{PlayerMajor}`
# holds and the `PayoffVector{L}(...)` constructors below apply to the alias.

"""
    GAMPayoffVector{N,T}

Alias for `PayoffVector{PlayerMajor,N,T}`: payoff values are ordered as in the
GameTracer .gam format:
1. Player-major blocks: player 1, ..., player N.
2. Within each block, action profiles are ordered with player 1 varying fastest,
   then player 2, ..., player N (i.e., column-major order).
"""
const GAMPayoffVector{N,T<:Real} = PayoffVector{PlayerMajor,N,T}

"""
    NFGPayoffVector{N,T}

Alias for `PayoffVector{ProfileMajor,N,T}`: payoff values are ordered as in the
Gambit .nfg format:
1. Profile-major blocks: action profiles are ordered with player 1 varying
   fastest, then player 2, ..., player N (i.e., column-major order).
2. Within each block, the payoffs to player 1, ..., player N.
"""
const NFGPayoffVector{N,T<:Real} = PayoffVector{ProfileMajor,N,T}

num_players(::PayoffVector{L,N}) where {L,N} = N

PayoffVector{L}(
    nums_actions::NTuple{N,Int}, payoffs::Vector{T}
) where {L<:PayoffLayout,N,T<:Real} = PayoffVector{L,N,T}(nums_actions, payoffs)

PayoffVector{L}(
    ::Type{T}, nums_actions::NTuple{N,Int}, payoffs::AbstractVector
) where {L<:PayoffLayout,N,T<:Real} =
    PayoffVector{L,N,T}(nums_actions, convert(Vector{T}, payoffs))
PayoffVector{L}(
    nums_actions::NTuple{N,Int}, payoffs::AbstractVector{T}
) where {L<:PayoffLayout,N,T<:Real} = PayoffVector{L}(T, nums_actions, payoffs)


# Forward: (i, i+1, ..., N, 1, ..., i-1)
@inline _perm_fwd(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(i + k - 1, N), Val(N))

# Backward: inverse rotation
@inline _perm_back(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(k - i + 1, N), Val(N))

@inline _colons(::Val{N}) where {N} = ntuple(_ -> Colon(), Val(N))

"""
    _player_block(p, i)

Return a view of `p.payoffs` holding the payoffs to player `i`, as an `N`-dim
array indexed by the action profile `(a_1, ..., a_N)`. This is the only place
where the layout `L` matters; no copy is made.
"""
@inline _player_block(p::PayoffVector{PlayerMajor,N}, i::Int) where {N} =
    view(reshape(p.payoffs, (p.nums_actions..., N)), _colons(Val(N))..., i)
@inline _player_block(p::PayoffVector{ProfileMajor,N}, i::Int) where {N} =
    view(reshape(p.payoffs, (N, p.nums_actions...)), i, _colons(Val(N))...)


"""
    PayoffVector{L}([T], g)

Construct a `PayoffVector` of layout `L` (and of eltype `T` if specified) from
a NormalFormGame `g`. `GAMPayoffVector([T], g)` and `NFGPayoffVector([T], g)`
are the versions for the two layouts.

# Examples

```julia
julia> player1 = Player([1 4; 2 5; 3 6]);

julia> player2 = Player([7 8 9; 10 11 12]);

julia> g = NormalFormGame(player1, player2)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)

julia> p = GameTheory.GAMPayoffVector(g);

julia> @show p.payoffs;
p.payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

julia> p = GameTheory.NFGPayoffVector(g);

julia> @show p.payoffs;
p.payoffs = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]
```
"""
function PayoffVector{L}(
    ::Type{T}, g::NormalFormGame{N}
) where {L<:PayoffLayout,N,T<:Real}
    nums_actions = g.nums_actions
    payoffs = Vector{T}(undef, prod(nums_actions)*N)
    p = PayoffVector{L,N,T}(nums_actions, payoffs)

    ntuple(Val(N)) do i
        copyto!(
            _player_block(p, i),
            PermutedDimsArray(g.players[i].payoff_array, _perm_back(Val(N), Val(i)))
        )
        nothing
    end

    return p
end

PayoffVector{L}(g::NormalFormGame{N,T}) where {L<:PayoffLayout,N,T<:Real} =
    PayoffVector{L}(T, g)


"""
    PayoffVector{L}([T], p)

Construct a `PayoffVector` of layout `L` (and of eltype `T` if specified) from
a PayoffVector `p` of possibly another layout. The payoffs are copied.

# Examples

```julia
julia> p = GameTheory.GAMPayoffVector((3, 2), collect(1:12));

julia> p_nfg = GameTheory.NFGPayoffVector(p);

julia> @show p_nfg.payoffs;
p_nfg.payoffs = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]
```
"""
function PayoffVector{L}(
    ::Type{T}, p::PayoffVector{L1,N}
) where {L<:PayoffLayout,L1<:PayoffLayout,N,T<:Real}
    payoffs = Vector{T}(undef, length(p.payoffs))
    p_new = PayoffVector{L,N,T}(p.nums_actions, payoffs)

    ntuple(Val(N)) do i
        copyto!(_player_block(p_new, i), _player_block(p, i))
        nothing
    end

    return p_new
end

PayoffVector{L}(p::PayoffVector{L1,N,T}) where {L<:PayoffLayout,L1<:PayoffLayout,N,T<:Real} =
    PayoffVector{L}(T, p)


"""
    NormalFormGame([T], p)

Construct a NormalFormGame (of eltype `T` if specified) from a PayoffVector
`p`.

# Examples

```julia
julia> nums_actions = (3, 2);

julia> payoffs = collect(1:12);

julia> @show payoffs;
payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

julia> p = GameTheory.GAMPayoffVector(nums_actions, payoffs);

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)

julia> p = GameTheory.NFGPayoffVector(nums_actions, payoffs);

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 2)  (7, 8)
 (3, 4)  (9, 10)
 (5, 6)  (11, 12)
```
"""
function NormalFormGame(::Type{T}, p::PayoffVector{L,N}) where {L<:PayoffLayout,N,T<:Real}
    players = ntuple(Val(N)) do i
        Player(
            T,
            PermutedDimsArray(_player_block(p, i), _perm_fwd(Val(N), Val(i)))
        )
    end

    return NormalFormGame{N,T}(players, p.nums_actions)
end

NormalFormGame(p::PayoffVector{L,N,T}) where {L<:PayoffLayout,N,T<:Real} =
    NormalFormGame(T, p)
