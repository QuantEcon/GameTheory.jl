#=
Compute all mixed Nash equilibria of a 2-player (non-degenerate) normal
form game by support enumeration.

Julia version of QuantEcon.py/support_enumeration.py

Authors: Daisuke Oyama, Zejin Shi

References
----------
B. von Stengel, "Equilibrium Computation for Two-Player Games in
Strategic and Extensive Form," Chapter 3, N. Nisan, T. Roughgarden, E.
Tardos, and V. Vazirani eds., Algorithmic Game Theory, 2007.
=#
import Compat.LinearAlgebra: LAPACKException, SingularException
import QuantEcon: next_k_array!

# For v0.6 compatibility
@static if !isdefined(Compat.LinearAlgebra, :lu!)
    lu!(A::StridedMatrix) = lufact!(A)
end

@static if !isdefined(Compat.LinearAlgebra, :ldiv!)
    ldiv!(A::LU{<:Any, <:StridedMatrix}, B::StridedVecOrMat) = A_ldiv_B!(A, B)
end

"""
    support_enumeration(g)

Compute mixed-action Nash equilibria with equal support size
for a 2-player normal form game by support enumeration. For a
non-degenerate game input, these are all the Nash equilibria.

The algorithm checks all the equal-size support pairs; if the
players have the same number n of actions, there are 2n choose n
minus 1 such pairs. This should thus be used only for small games.

# Arguments

- `g::NormalFormGame{2,T}`: 2-player NormalFormGame instance.

# Returns

- `::Vector{NTuple{2,Vector{S}}}`: Mixed-action Nash equilibria that are found,
  where `S` is Float if `T` is Int or Float, and Rational if `T` is Rational.
"""
function support_enumeration(g::NormalFormGame{2,T}) where T
    S = typeof(zero(T)/one(T))
    c = Channel{Tuple{Vector{S},Vector{S}}}(0)
    task = support_enumeration_task(c, g)
    bind(c, task)
    schedule(task)
    NEs = collect(c)

    return NEs

end

"""
    support_enumeration_task(c, g)

Task version of `support_enumeration`.

# Arguments

- `c::Channel`: Channel to be binded with the support enumeration task.
- `g::NormalFormGame{2}`: 2-player NormalFormGame instance.

# Returns

- `::Task`: Runnable task for generating Nash equilibria.
"""
function support_enumeration_task(c::Channel,
                                  g::NormalFormGame{2})

    task = Task(
        () -> _support_enumeration_producer(c,
                                            (g.players[1].payoff_array,
                                             g.players[2].payoff_array))
    )

    return task
end

"""
    _support_enumeration_producer(c, payoff_matrices)

Main body of `support_enumeration_task`.

# Arguments

- `c::Channel`: Channel to be binded with the support enumeration task.
- `payoff_matrices::NTuple{2,Matrix{T}}`: Payoff matrices of player 1 and
  player 2, where `T<:Real`.

# Puts

- `NTuple{2,Vector{S}}`: Tuple of Nash equilibrium mixed actions, where `S` is
  Float if `T` is Int or Float, and Rational if `T` is Rational.
"""
function _support_enumeration_producer(c::Channel,
                                       payoff_matrices
                                       ::NTuple{2,Matrix{T}}) where T<:Real

    nums_actions = size(payoff_matrices[1], 1), size(payoff_matrices[2], 1)
    n_min = min(nums_actions...)
    S = typeof(zero(T)/one(T))

    for k = 1:n_min
        supps = (collect(1:k), Vector{Int}(undef, k))
        actions = (Vector{S}(undef, k), Vector{S}(undef, k))
        A = Matrix{S}(undef, k+1, k+1)
        b = Vector{S}(undef, k+1)
        while supps[1][end] <= nums_actions[1]
            supps[2][:] = collect(1:k)
            while supps[2][end] <= nums_actions[2]
                if _indiff_mixed_action!(A, b, actions[2],
                                         payoff_matrices[1],
                                         supps[1], supps[2])
                    if _indiff_mixed_action!(A, b, actions[1],
                                             payoff_matrices[2],
                                             supps[2], supps[1])
                        out = (zeros(S, nums_actions[1]),
                               zeros(S, nums_actions[2]))
                        for (p, (supp, action)) in enumerate(zip(supps,
                                                                 actions))
                            out[p][supp] = action
                        end
                        put!(c, out)
                    end
                end
                next_k_array!(supps[2])
            end
            next_k_array!(supps[1])
        end
    end

end

function _solve!(A::Matrix{T}, b::Vector{T}) where T <: Union{Float64,Float32}
    r = 0
    try
        LAPACK.gesv!(A, b)
    catch LAPACKException
        r = 1
    end
    return r
end

@inline function _solve!(A::Matrix{Rational{T}},
                         b::Vector{Rational{T}}) where T <: Integer
    r = 0
    try
        b[:] = ldiv!(lu!(A), b)
    catch SingularException
        r = 1
    end
    return r
end

"""
    _indiff_mixed_action!(A, b, out, payoff_matrix, own_supp, opp_supp)

Given a player's payoff matrix `payoff_matrix`, an array `own_supp`
of this player's actions, and an array `opp_supp` of the opponent's
actions, each of length k, compute the opponent's mixed action whose
support equals `opp_supp` and for which the player is indifferent
among the actions in `own_supp`, if any such exists. Return `true`
if such a mixed action exists and actions in `own_supp` are indeed
best responses to it, in which case the outcome is stored in `out`;
`false` otherwise. Arrays `A` and `b` are used in intermediate
steps.

# Arguments

- `A::Matrix{T}`: Matrix of shape (k+1, k+1) used in intermediate steps, where
  `T<:Real`.
- `b::Vector{T}`: Vector of length k+1 used in intermediate steps, where
  `T<:Real`.
- `out::Vector{T}`: Vector of length k to store the nonzero values of the
  desired mixed action, where `T<:Real`.
- `payoff_matrix::Matrix`: The player's payoff matrix, of shape (m, n).
- `own_supp::Vector{Int}`: Vector containing the player's action indices, of
  length k.
- `opp_supp::Vector{Int}`: Vector containing the opponent's action indices, of
  length k.

# Returns

- `::Bool`: `true` if a desired mixed action exists and `false` otherwise.
"""
function _indiff_mixed_action!(A::Matrix{T}, b::Vector{T},
                               out::Vector{T},
                               payoff_matrix::Matrix,
                               own_supp::Vector{Int},
                               opp_supp::Vector{Int}) where T<:Real

    m = size(payoff_matrix, 1)
    k = length(own_supp)

    for j in 1:k, i in 1:k
        A[i, j] = payoff_matrix[own_supp[i], opp_supp[j]]
    end
    A[1:end-1, end] .= -one(T)
    A[end, 1:end-1] .= one(T)
    A[end, end] = zero(T)
    b[1:end-1] .= zero(T)
    b[end] = one(T)

    r = _solve!(A, b)
    r == 0 || return false  # A: singular

    for i in 1:k
        b[i] <= zero(T) && return false
    end

    out[:] = b[1:end-1]
    val = b[end]

    if k == m
        return true
    end

    own_supp_flags = falses(m)
    own_supp_flags[own_supp] .= true

    for i = 1:m
        if !own_supp_flags[i]
            payoff = zero(T)
            for j = 1:k
                payoff += payoff_matrix[i, opp_supp[j]] * out[j]
            end
            if payoff > val
                return false
            end
        end
    end

    return true
end
