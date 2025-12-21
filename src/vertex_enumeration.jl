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


struct BestResponsePolytope{T,PT<:Polyhedron{T}}
    ndim::Int
    poly::PT
end

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

    return BestResponsePolytope{S,typeof(poly)}(m, poly)
end


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


function _vertex_enumeration_producer(
        c::Channel, brps::NTuple{2,BestResponsePolytope{T,PT}}
    ) where {T,PT}
    m, n = brps[1].ndim, brps[2].ndim
    ZERO_LABELING1_BITS = (UInt64(1) << UInt64(m)) - UInt64(1)
    COMPLETE_LABELING_BITS = (UInt64(1) << UInt64(m + n)) - UInt64(1)

    vertices2 = points(vrep(brps[2].poly))
    labelings_bits_dict2 = Dict{UInt64,Vector{T}}()
    sizehint!(labelings_bits_dict2, length(vertices2))
    for (idx, pt) in zip(eachindex(vertices2), vertices2)
        indices = incidenthalfspaceindices(brps[2].poly, idx)
        labelings_bits_dict2[_indices_to_bits(indices)] = pt
    end

    vertices1 = points(vrep(brps[1].poly))
    for (idx, pt1) in zip(eachindex(vertices1), vertices1)
        indices = incidenthalfspaceindices(brps[1].poly, idx)
        bits1 = _indices_to_bits(indices)
        bits1 == ZERO_LABELING1_BITS && continue

        complement1 = xor(bits1, COMPLETE_LABELING_BITS)
        pt2 = get(labelings_bits_dict2, complement1, nothing)
        pt2 === nothing && continue
        put!(c, (pt1/sum(pt1), pt2/sum(pt2)))
    end
end


"""
    _indices_to_bits(indices)

Convert a vector of `Polyhedra.AbstractIndex` representing the set bits into
the corresponding integer.

# Arguments

- `indices::Vector{<:Polyhedra.AbstractIndex}`: Vector of `Polyhedra.AbstractIndex`
  with distinct integer values from 1, ..., 64.

# Returns

- `bits::UInt64`: Integer with set bits represented by the input integers.

# Examples

```julia
julia> indices = map(Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}, [1, 2, 3])
3-element Vector{Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}}:
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(1)
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(2)
 Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(3)

julia> _indices_to_bits(indices)
0x0000000000000007
```
"""
function _indices_to_bits(indices::Vector{<:Polyhedra.AbstractIndex})
    bits = UInt64(0)
    @inbounds for idx in indices
        bits |= UInt64(1) << UInt64(idx.value-1)
    end
    return bits
end
