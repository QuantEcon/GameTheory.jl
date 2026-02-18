# GamPayoffVector #

"""
    GamPayoffVector{N,T}

Internal intermediate representation that stores payoffs in a single flat
vector.

Payoff values are ordered as in the GameTracer .gam format:
1. Player-major blocks: player 1, ..., player N.
2. Within each block, action profiles are ordered with player 1 varying fastest,
   then player 2, ..., player N (i.e., Fortran/column-major order).

# Fields

- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `payoffs::Vector{T}` : Vector storing payoffs in .gam order.
"""
struct GamPayoffVector{N,T<:Real}
    nums_actions::NTuple{N,Int}
    payoffs::Vector{T}

    function GamPayoffVector{N,T}(
        nums_actions::NTuple{N,Int}, payoffs::Vector{T}
    ) where {N,T<:Real}
        any(n -> n <= 0, nums_actions) &&
            throw(ArgumentError("all nums_actions must be positive"))
        expected = prod(nums_actions) * N
        length(payoffs) == expected || throw(ArgumentError(
            "payoffs length mismatch: expected $expected, got $(length(payoffs))"
        ))
        return new(nums_actions, payoffs)
    end
end

num_players(::GamPayoffVector{N}) where {N} = N

GamPayoffVector(
    nums_actions::NTuple{N,Int}, payoffs::Vector{T}
) where {N,T<:Real} = GamPayoffVector{N,T}(nums_actions, payoffs)
GamPayoffVector(
    nums_actions::NTuple{N,Int}, payoffs::AbstractVector{T}
) where {N,T<:Real} =
    GamPayoffVector{N,T}(nums_actions, convert(Vector{T}, payoffs))


# Forward: (i, i+1, ..., N, 1, ..., i-1)
@inline perm_fwd(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(i + k - 1, N), Val(N))

# Backward: inverse rotation
@inline perm_back(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(k - i + 1, N), Val(N))


"""
    GamPayoffVector([T], g)

Construct a GamPayoffVector (of eltype `T` if specified) from a NormalFormGame
`g`.

# Examples

```julia
julia> player1 = Player([1 4; 2 5; 3 6]);

julia> player2 = Player([7 8 9; 10 11 12]);

julia> g = NormalFormGame(player1, player2);

julia> println(g)
3×2 NormalFormGame{2, Int64}:
 [1, 7]  [4, 10]
 [2, 8]  [5, 11]
 [3, 9]  [6, 12]

julia> p = GamPayoffVector(g);

julia> @show p.payoffs;
p.payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
```
"""
function GamPayoffVector(::Type{T}, g::NormalFormGame{N}) where {N,T<:Real}
    nums_actions = g.nums_actions
    na = prod(nums_actions)
    payoffs = Vector{T}(undef, na*N)

    ntuple(Val(N)) do i
        copyto!(
            reshape(view(payoffs, na*(i-1)+1:na*i), nums_actions),
            PermutedDimsArray(g.players[i].payoff_array, perm_back(Val(N), Val(i)))
        )
        nothing
    end

    return GamPayoffVector{N,T}(nums_actions, payoffs)
end

GamPayoffVector(g::NormalFormGame{N,T}) where {N,T<:Real} = GamPayoffVector(T, g)


"""
    NormalFormGame([T], p)

Construct a NormalFormGame (of eltype `T` if specified) from a GamPayoffVector
`p`.

# Examples

```julia
julia> nums_actions = (3, 2);

julia> payoffs = collect(1:12);

julia> @show payoffs;
payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

julia> p = GamPayoffVector(nums_actions, payoffs);

julia> g = NormalFormGame(p);

julia> println(g)
3×2 NormalFormGame{2, Int64}:
 [1, 7]  [4, 10]
 [2, 8]  [5, 11]
 [3, 9]  [6, 12]
```
"""
function NormalFormGame(::Type{T}, p::GamPayoffVector{N}) where {N,T<:Real}
    nums_actions = p.nums_actions
    na = prod(nums_actions)

    players = ntuple(Val(N)) do i
        Player(
            T,
            PermutedDimsArray(
                reshape(view(p.payoffs, na*(i-1)+1:na*i), nums_actions),
                perm_fwd(Val(N), Val(i))
            )
        )
    end

    return NormalFormGame{N,T}(players, nums_actions)
end

NormalFormGame(p::GamPayoffVector{N,T}) where {N,T<:Real} = NormalFormGame(T, p)
