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

"""
    support_enumeration(g)

Compute mixed-action Nash equilibria with equal support size
for a 2-player normal form game by support enumeration. For a
non-degenerate game input, these are all the Nash equilibria.

The algorithm checks all the equal-size support pairs; if the
players have the same number n of actions, there are 2n choose n
minus 1 such pairs. This should thus be used only for small games.

# Arguments

- `g::NormalFormGame{2}`: 2-player NormalFormGame instance.

# Returns

- `::Vector{Tuple{Vector{Real}, Vector{Real}}}`: Mixed-action
  Nash equilibria that are found.
"""
function support_enumeration(g::NormalFormGame{2})

    c = Channel(0)
    task = support_enumeration_task(c, g)
    bind(c, task)
    schedule(task)
    NEs = Tuple{Vector{Real}, Vector{Real}}[NE for NE in c]

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
- `payoff_matrices::NTuple{2, Matrix{T}}`: Payoff matrices of player 1 and
  player 2. T<:Real.

# Puts

- `Tuple{Vector{S},Vector{S}}`: Tuple of Nash equilibrium mixed actions.
  `S` is Float if `T` is Int or Float, and Rational if `T` is Rational.
"""
function _support_enumeration_producer{T<:Real}(c::Channel,
                                                payoff_matrices
                                                ::NTuple{2,Matrix{T}})

    nums_actions = size(payoff_matrices[1], 1), size(payoff_matrices[2], 1)
    n_min = min(nums_actions...)
    S = typeof(zero(T)/one(T))

    for k = 1:n_min
        supps = (collect(1:k), Vector{Int}(k))
        actions = (Vector{S}(k), Vector{S}(k))
        A = Matrix{S}(k+1, k+1)
        b = Vector{S}(k+1)
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
                _next_k_array!(supps[2])
            end
            _next_k_array!(supps[1])
        end
    end

end

"""
    _indiff_mixed_action!(A, b, out, payoff_matrix,
                          own_supp, opp_supp)

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

- `A::Matrix{T}`: Matrix used in intermediate steps. T<:Real.
- `b::Vector{T}`: Vector used in intermediate steps. T<:Real.
- `out::Vector{T}`: Vector to store the nonzero values of the
  desired mixed action. T<:Real.
- `payoff_matrix::Matrix`: The player's payoff matrix.
- `own_supp::Vector{Int}`: Vector containing the player's action indices.
- `opp_supp::Vector{Int}`: Vector containing the opponent's action indices.

# Returns

- `::Bool`: `true` if a desired mixed action exists and `false` otherwise.
"""
function _indiff_mixed_action!{T<:Real}(A::Matrix{T}, b::Vector{T},
                                        out::Vector{T},
                                        payoff_matrix::Matrix,
                                        own_supp::Vector{Int},
                                        opp_supp::Vector{Int})

    m = size(payoff_matrix, 1)
    k = length(own_supp)

    A[1:end-1, 1:end-1] = payoff_matrix[own_supp, opp_supp]
    A[1:end-1, end] = -one(T)
    A[end, 1:end-1] = one(T)
    A[end, end] = zero(T)
    b[1:end-1] = zero(T)
    b[end] = one(T)
    try
        b = A_ldiv_B!(lufact!(A), b)
    catch LinAlg.SingularException
        return false
    end

    for i in 1:k
        b[i] <= zero(T) && return false
    end

    out[:] = b[1:end-1]
    val = b[end]

    if k == m
        return true
    end

    own_supp_flags = falses(m)
    own_supp_flags[own_supp] = true

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

"""
    _next_k_combination(x)

Find the next k-combination, as described by an integer in binary
representation with the k set bits, by "Gosper's hack".

Copy-paste from en.wikipedia.org/wiki/Combinatorial_number_system

# Arguments

- `x::Int`: Integer with k set bits.

# Returns

- `::Int`: Smallest integer > x with k set bits.
"""
function _next_k_combination(x::Int)

    u = x & -x
    v = u + x
    return v + (fld((v ⊻ x), u) >> 2)

end

"""
    _next_k_array!(a)

Given an array `a` of k distinct nonnegative integers, return the
next k-array in lexicographic ordering of the descending sequences
of the elements. `a` is modified in place.

# Arguments

- `a::Vector{Int}`: Array of length k.

# Returns

- `:::Vector{Int}`: Next k-array of `a`.

# Examples

```julia
julia> n, k = 4, 2
(4,2)

julia> a = collect(1:k)
2-element Array{Int64,1}:
 1
 2

julia> while a[end] < n + 1
           @show a
           _next_k_array!(a)
       end
a = [1,2]
a = [1,3]
a = [2,3]
a = [1,4]
a = [2,4]
a = [3,4]
```
"""
function _next_k_array!(a::Vector{Int})

    k = length(a)
    if k == 0
        return a
    end

    x = 0
    for i = 1:k
        x += (1 << (a[i] - 1))
    end

    x = _next_k_combination(x)

    pos = 0
    for i = 1:k
        while x & 1 == 0
            x = x >> 1
            pos += 1
        end
        a[i] = pos + 1
        x = x >> 1
        pos += 1
    end

    return a
end
