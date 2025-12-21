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

struct BestResponsePolytope{T,PT<:Polyhedron{T}}
    poly::PT
end

_coeftype(::Type{<:Rational}) = Rational{BigInt}
_coeftype(::Type{<:AbstractFloat}) = Float64
_coeftype(::Type{<:Integer}) = Float64

function BestResponsePolytope(opponent_player::Player{2,T}; idx=1,
                              plib::Polyhedra.Library=CDDLib.Library()) where T
    S = typeof(zero(T)/one(T))
    S2 = _coeftype(S)
    B = opponent_player.payoff_array
    n, m = size(B)
    lib = Polyhedra.similar_library(plib, m, S2)
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
    hrep(poly)

    return BestResponsePolytope{S2,typeof(poly)}(poly)
end
