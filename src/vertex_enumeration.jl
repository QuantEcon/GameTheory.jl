#=
Compute all mixed Nash equilibria of a 2-player normal form game by
vertex enumeration.

Julia version of QuantEcon.py/vertex_enumeration.py

References
----------
B. von Stengel, "Equilibrium Computation for Two-Player Games in
Strategic and Extensive Form," Chapter 3, N. Nisan, T. Roughgarden, E.
Tardos, and V. Vazirani eds., Algorithmic Game Theory, 2007.
=#

using LinearAlgebra: diagind

_coeftype(::Type{<:Rational}) = Rational{BigInt}
_coeftype(::Type{<:AbstractFloat}) = Float64
_coeftype(::Type{<:Integer}) = Float64


@doc doc"""
    BestResponsePolytope{T,PT}

Type that represents the best response polytope for a player in a two-player
normal form game.

Let ``A`` and ``B`` be the m x n and n x m payoff matrices of players 1 and 2,
respectively, where the payoffs are assumed to have been shifted in such a way
that ``A`` and ``B`` are nonnegative and have no zero column. In von Stengel
(2007), the best response polytope for player 1 is defined by

```math
P = \{x \in \mathbb{R}^m \mid x \geq 0,\ B x \leq 1\},
```

and that for player 2 by

```math
Q = \{y \in \mathbb{R}^n \mid A y \leq 1,\ y \geq 0\}.
```

The polytope is stored in the field `poly` as a `Polyhedra.Polyhedron{T}` which
conducts vertex enumeration.

# Fields

- `idx::Int`: Player index in the normal form game, either 1 or 2.
- `ndim::Int`: Dimension of the polytope.
- `poly::PT`: `Polyhedra.Polyhedron{T}` instance representing the best response
  polytope.
"""
struct BestResponsePolytope{T,PT<:Polyhedron{T}}
    idx::Int
    ndim::Int
    poly::PT
end


"""
    BestResponsePolytope(opponent_player; idx=1, plib=CDDLib.Library())

Constructor of a `BestResponsePolytope` for player `idx`.

# Arguments

- `opponent_player::Player{2,T}`: Instance of Player with one opponent.
- `idx`: Player index in the normal form game, either 1 or 2.
- `plib::Polyhedra.Library`: Allows to choose a particular package for
  polyhedral computation.

# Returns

- `::BestResponsePolytope{S,PT}`: `BestResponsePolytope` instance, where `S` is
  `Float64` if `T` is Int or Float, and `Rational{BigInt}` if `T` is Rational.
"""
function BestResponsePolytope(opponent_player::Player{2,T}; idx=1,
                              plib::Polyhedra.Library=CDDLib.Library()) where T
    S = _coeftype(T)
    B = opponent_player.payoff_array
    n, m = size(B)
    lib = Polyhedra.similar_library(plib, m, S)
    D = Matrix{S}(undef, m + n, m)

    nonneg_cond_start, payoff_cond_start = idx == 1 ? (1, m + 1) : (n + 1, 1)

    col_mins = vec(minimum(B, dims = 1))
    col_maxs = vec(maximum(B, dims = 1))
    nonpos_const_cols = (col_maxs .== col_mins) .& (col_mins .<= 0)
    shifts = zeros(T, m)
    negcols = col_mins .< 0
    shifts[negcols] .= -col_mins[negcols]
    shifts[nonpos_const_cols] .+= 1
    D_payoff = @view D[payoff_cond_start : payoff_cond_start + n - 1, :]
    @inbounds for j in 1:m
        sj = shifts[j]
        for i in 1:n
            D_payoff[i, j] = B[i, j] + sj
        end
    end

    rows = nonneg_cond_start : nonneg_cond_start + m - 1
    D[rows, :] .= 0
    D_nonneg = @view D[rows, :]
    D_nonneg[diagind(D_nonneg)] .= -1

    b = zeros(S, m + n)
    b[payoff_cond_start : payoff_cond_start + n - 1] .= 1

    hr = hrep(D, b)
    poly = polyhedron(hr, lib)

    return BestResponsePolytope{S,typeof(poly)}(idx, m, poly)
end


"""
    vertex_enumeration(g; plib=CDDLib.Library())

Compute mixed-action Nash equilibria of a 2-player normal form game by
enumeration and matching of vertices of the best response polytopes. For
a non-degenerate game input, these are all the Nash equilibria.

Internally, the library specified by `plib` is used to compute vertex
enumeration of the best response polytopes. Then, for each vertex of the
polytope for player 1, vertices of the polytope for player 2 are searched
to find a completely labeled pair.

# Arguments

- `g::NormalFormGame{2,T}`: 2-player NormalFormGame instance.
- `plib::Polyhedra.Library`: Allows to choose a particular package for
  polyhedral computation.

# Returns

- `::Vector{NTuple{2,Vector{S}}}`: Mixed-action Nash equilibria that are found,
  where `S` is `Float64` if `T` is Int or Float, and `Rational{BigInt}` if `T`
  is Rational.

# Examples

```julia
julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> player1 = Player([3 3; 2 5; 0 6]);

julia> player2 = Player([3 2 3; 2 6 1]);

julia> g = NormalFormGame(player1, player2);

julia> println(g)
3×2 NormalFormGame{2, Int64}:
 [3, 3]  [3, 2]
 [2, 2]  [5, 6]
 [0, 3]  [6, 1]

julia> vertex_enumeration(g)
3-element Vector{Tuple{Vector{Float64}, Vector{Float64}}}:
 ([1.0, 0.0, 0.0], [1.0, 0.0])
 ([0.8, 0.2, 0.0], [0.666667, 0.333333])
 ([0.0, 0.333333, 0.666667], [0.333333, 0.666667])
```
"""
function vertex_enumeration(g::NormalFormGame{2,T};
                            plib::Polyhedra.Library=CDDLib.Library()) where T
    S = _coeftype(T)
    c = Channel{Tuple{Vector{S},Vector{S}}}(0)
    task = vertex_enumeration_task(c, g, plib=plib)
    bind(c, task)
    schedule(task)
    NEs = collect(c)

    return NEs
end


"""
    vertex_enumeration_task(c, g; plib=CDDLib.Library())

Task version of `vertex_enumeration`.

# Arguments

- `c::Channel`: Channel to be binded with the support enumeration task.
- `g::NormalFormGame{2}`: 2-player NormalFormGame instance.
- `plib::Polyhedra.Library`: Allows to choose a particular package for
  polyhedral computation.

# Returns

- `::Task`: Runnable task for generating Nash equilibria.

# Examples

```julia
julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> player1 = Player([3 3; 2 5; 0 6]);

julia> player2 = Player([3 2 3; 2 6 1]);

julia> g = NormalFormGame(player1, player2);

julia> println(g)
3×2 NormalFormGame{2, Int64}:
 [3, 3]  [3, 2]
 [2, 2]  [5, 6]
 [0, 3]  [6, 1]

julia> c = Channel{Tuple{Vector{Float64},Vector{Float64}}}(0);

julia> t = vertex_enumeration_task(c, g);

julia> bind(c, t); schedule(t);

julia> for NE in c
           display(NE)
       end
([1.0, 0.0, 0.0], [1.0, 0.0])
([0.8, 0.2, 0.0], [0.666667, 0.333333])
([0.0, 0.333333, 0.666667], [0.333333, 0.666667])
```
"""
function vertex_enumeration_task(c::Channel,
                                 g::NormalFormGame{2};
                                 plib::Polyhedra.Library=CDDLib.Library())
    brps = ntuple(
        i -> BestResponsePolytope(g.players[3-i], idx=i, plib=plib), 2
    )
    task = Task(
        () -> _vertex_enumeration_producer(c, brps)
    )

    return task
end


mutable struct LabelBitsIterator{T,PT<:Polyhedron{T}}
    brp::BestResponsePolytope{T,PT}
    num_vertices::Int
    hincidence::Vector{BitSet}
    warned::Bool
end

_copyincidence(poly::CDDLib.Polyhedron) = copyincidence(CDDLib.getpoly(poly))

function LabelBitsIterator(
    brp::BestResponsePolytope{T,PT}
) where {T,PT<:CDDLib.Polyhedron{T}}
    hincidence = _copyincidence(brp.poly)
    num_vertices = length(hincidence)
    return LabelBitsIterator{T,PT}(brp, num_vertices, hincidence, false)
end

function LabelBitsIterator(brp::BestResponsePolytope{T,PT}) where {T,PT}
    num_vertices = npoints(vrep(brp.poly))
    hincidence = BitSet[]
    return LabelBitsIterator{T,PT}(brp, num_vertices, hincidence, false)
end

Base.IteratorEltype(::Type{<:LabelBitsIterator}) = Base.HasEltype()
Base.eltype(::Type{<:LabelBitsIterator}) = UInt64

Base.IteratorSize(::Type{<:LabelBitsIterator}) = Base.HasLength()
Base.length(it::LabelBitsIterator) = it.num_vertices

function Base.iterate(
    it::LabelBitsIterator{T,<:CDDLib.Polyhedron{T}}, state::Int=1
) where {T}
    state > it.num_vertices && return nothing

    indices = it.hincidence[state]
    num_indices = length(indices)
    ndim = it.brp.ndim
    if num_indices > ndim
        if it.brp.idx == 1
            _delete_largest_k!(indices, num_indices-ndim)
        else
            _delete_smallest_k!(indices, num_indices-ndim)
        end
        if !it.warned
            @warn "Payoff degeneracy detected for Player $(it.brp.idx)"
            it.warned = true
        end
    end
    bits = indices.bits[1] >> 1
    return (bits, state+1)
end

function Base.iterate(it::LabelBitsIterator{T,PT}, state::Int=1) where {T,PT}
    state > it.num_vertices && return nothing

    idx = Polyhedra.Index{T,Vector{T}}(state)
    indices = incidenthalfspaceindices(it.brp.poly, idx)
    num_indices = length(indices)
    ndim = it.brp.ndim
    if num_indices > ndim
        if it.brp.idx == 1
            bits = _indices_to_bits(@view indices[1:ndim])
        else
            bits = _indices_to_bits(@view indices[end-ndim+1:end])
        end
        if !it.warned
            @warn "Payoff degeneracy detected for Player $(it.brp.idx)"
            it.warned = true
        end
    else
        bits = _indices_to_bits(indices)
    end
    return (bits, state+1)
end


"""
    _vertex_enumeration_producer(c, brps)

Main body of `vertex_enumeration_task`.

# Arguments

- `c::Channel`: Channel to be binded with the support enumeration task.
- `brps::NTuple{2,BestResponsePolytope{T,PT}}`: Tuple of Best response
  polytopes of player 1 and player 2, where `PT<:Polyhedron{T}`.

# Puts

- `NTuple{2,Vector{T}}`: Tuple of Nash equilibrium mixed actions.
"""
function _vertex_enumeration_producer(
        c::Channel, brps::NTuple{2,BestResponsePolytope{T,PT}}
    ) where {T,PT}
    m, n = brps[1].ndim, brps[2].ndim
    @assert m + n <= 63 "Implemented only for games with sum(g.nums_actions) <= 63"

    ZERO_LABELING2_BITS = ((UInt64(1) << UInt64(n)) - UInt64(1)) << UInt64(m)
    COMPLETE_LABELING_BITS = (UInt64(1) << UInt64(m + n)) - UInt64(1)

    it2 = LabelBitsIterator(brps[2])
    labelings_bits_dict2 = Dict{UInt64,Int}()
    sizehint!(labelings_bits_dict2, length(it2))
    @inbounds for (j, bits2) in enumerate(it2)
        labelings_bits_dict2[bits2] = j
    end
    delete!(labelings_bits_dict2, ZERO_LABELING2_BITS)

    it1 = LabelBitsIterator(brps[1])
    @inbounds for (i, bits1) in enumerate(it1)
        complement1 = xor(bits1, COMPLETE_LABELING_BITS)
        j = get(labelings_bits_dict2, complement1, 0)
        j == 0 && continue
        pt1 = get(brps[1].poly, Polyhedra.Index{T,Vector{T}}(i))
        pt2 = get(brps[2].poly, Polyhedra.Index{T,Vector{T}}(j))
        put!(c, (pt1/sum(pt1), pt2/sum(pt2)))
    end
end


function _clear_lowest_k_setbits(x::UInt64, k::Integer)
    for _ in 1:k
        x &= x - 1
        x == 0 && break
    end
    return x
end

function _delete_smallest_k!(indices::BitSet, k::Integer)
    indices.bits[1] = _clear_lowest_k_setbits(indices.bits[1], k)
    return indices
end

function _delete_largest_k!(indices::BitSet, k::Integer)
    indices.bits[1] = bitreverse(
        _clear_lowest_k_setbits(bitreverse(indices.bits[1]), k)
    )
    return indices
end


"""
    _indices_to_bits(indices)

Convert a vector of `Polyhedra.AbstractIndex` representing the set bits into
the corresponding integer.

# Arguments

- `indices::AbstractVector{<:Polyhedra.AbstractIndex}`: Vector of `Polyhedra.AbstractIndex`
  with distinct integer values from 1, ..., 64.

# Returns

- `bits::UInt64`: Integer with set bits represented by the input integers.

# Examples

```julia
julia> indices = Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}.([1, 2, 3])
3-element Vector{Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}}:
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(1)
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(2)
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(3)

julia> _indices_to_bits(indices)
0x0000000000000007
```
"""
function _indices_to_bits(indices::AbstractVector{<:Polyhedra.AbstractIndex})
    bits = UInt64(0)
    @inbounds for idx in indices
        bits |= UInt64(1) << UInt64(idx.value-1)
    end
    return bits
end
