#=
Compute all mixed Nash equilibria of a 2-player (non-degenerate) normal
form game by support enumeration.

References
----------
B. von Stengel, "Equilibrium Computation for Two-Player Games in
Strategic and Extensive Form," Chapter 3, N. Nisan, T. Roughgarden, E.
Tardos, and V. Vazirani eds., Algorithmic Game Theory, 2007.
=#

"""
    support_enumeration(g::NormalFormGame)

Compute mixed-action Nash equilibria with equal support size for a
2-player normal form game by support enumeration. For a
non-degenerate game input, these are all Nash equilibria.

The algorithm checks all the equal-size support pairs; if the
players have the same number n of actions, there are 2n choose n
minus 1 such pairs. This should thus be used only for small games.

# Arguments
* `g::NormalFormGame`: NormalFormGame instance.

# Returns
* `::Vector{Tuple{Vector{Float64},Vector{Float64}}}`: Mixed-action
    Nash equilibria that are found.
"""
function support_enumeration(g::NormalFormGame)

    task = support_enumeration_gen(g)

    NEs = Tuple{Vector{Float64},Vector{Float64}}[NE for NE in task]

    return NEs

end

"""
    support_enumeration_gen(g::NormalFormGame)

Generator version of `support_enumeration`.

# Arguments
* `g::NormalFormGame`: NormalFormGame instance.

# Returns
* `::Task`: runnable task for generating Nash equilibria.
"""
function support_enumeration_gen(g::NormalFormGame)

    N = length(g.nums_actions)
    if N != 2
        throw(ArgumentError("Implemented only for 2-player games"))
    end

    task = Task(() -> _support_enumeration_gen(g.players[1].payoff_array,
                                               g.players[2].payoff_array))

    return task
end

"""
    _support_enumeration_gen{T<:Real}(payoff_matrix1::Matrix{T},
                                      payoff_matrix2::Matrix{T})

Main body of `support_enumeration_gen`.

# Arguments
* `payoff_matrix1::Matrix{T}`: Payoff matrix of player 1.
* `payoff_matrix2::Matrix{T}`: Payoff matrix of player 2.

# Produces
* `Tuple{Vector{Float64},Vector{Float64}}`: Tuple of Nash equilibrium
    mixed actions.
"""
function _support_enumeration_gen{T<:Real}(payoff_matrix1::Matrix{T},
                                           payoff_matrix2::Matrix{T})

    nums_actions = size(payoff_matrix1, 1), size(payoff_matrix2, 1)
    n_min = min(nums_actions...)

    for k = 1:n_min
        supps = (collect(1:k), Vector{Int}(k))
        actions = (Vector{Float64}(k), Vector{Float64}(k))
        A = Matrix{T}(k+1, k+1)
        A[1:end-1, end] = -one(T)
        A[end, 1:end-1] = one(T)
        A[end, end] = zero(T)
        b = zeros(T, k+1)
        b[end] = one(T)
        while supps[1][end] < nums_actions[1]+1
            supps[2][:] = collect(1:k)
            while supps[2][end] < nums_actions[2]+1
                if _indiff_mixed_action!(A, actions[2], b,
                                         payoff_matrix1, supps[1], supps[2])
                    if _indiff_mixed_action!(A, actions[1], b, payoff_matrix2,
                                             supps[2], supps[1])
                        out = (zeros(nums_actions[1]),
                               zeros(nums_actions[2]))
                        for (p, (supp, action)) in enumerate(zip(supps,
                                                                 actions))
                            out[p][supp] = action
                        end
                        produce(out)
                    end
                end
                _next_k_array!(supps[2])
            end
            _next_k_array!(supps[1])
        end
    end

end

"""
    _indiff_mixed_action!{T<:Real}(A::Matrix{T}, out::Vector{Float64},
                                   b::Vector{T},
                                   payoff_matrix::Matrix{T},
                                   own_supp::Vector{Int},
                                   opp_supp::Vector{Int})

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
* `A::Matrix{T}`: Matrix used in intermediate steps.
* `out::Vector{Float64}`: Vector to store the nonzero values of the 
    desired mixed action.
* `b::Vector{T}`: Vector used in intermediate steps.
* `payoff_matrix::Matrix{T}`: The player's payoff matrix.
* `own_supp::Vector{Int}`: Vector containing the player's action indices.
* `opp_supp::Vector{Int}`: Vector containing the opponent's action indices.

# Returns
* `::Bool`: `true` if a desired mixed action exists and `false` otherwise.
"""
function _indiff_mixed_action!{T<:Real}(A::Matrix{T}, out::Vector{Float64},
                                        b::Vector{T},
                                        payoff_matrix::Matrix{T},
                                        own_supp::Vector{Int},
                                        opp_supp::Vector{Int})

    m = size(payoff_matrix, 1)
    k = length(own_supp)

    A[1:end-1, 1:end-1] = payoff_matrix[own_supp, :][:, opp_supp]
    sol = Vector{Float64}(k+1)
    try
        sol[:] = A \ b
    catch LinAlg.SingularException
        return false
    end

    if any(sol[1:end-1] .<= 0)
        return false
    end
    out[:] = sol[1:end-1]
    val = sol[end]

    if k == m
        return true
    end

    own_supp_flags = falses(m)
    own_supp_flags[own_supp] = true

    for i = 1:m
        if !own_supp_flags[i]
            payoff = 0.
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
    _next_k_combination(x::Int)

Find the next k-combination, as described by an integer in binary
representation with the k set bits, by "Gosper's hack".

Copy-paste from en.wikipedia.org/wiki/Combinatorial_number_system

# Arguments
* `x::Int`: Integer with k set bits.

# Returns
* `::Int`: Smallest integer > x with k set bits.
"""
function _next_k_combination(x::Int)

    u = x & -x
    v = u + x
    return v + (fld((v $ x), u) >> 2)

end

"""
    _next_k_array!(a::Vector{Int})

Given an array `a` of k distinct nonnegative integers, return the
next k-array in lexicographic ordering of the descending sequences
of the elements. `a` is modified in place.

# Arguments
* `a::Vector{Int}`: Array of length k.

# Returns
* `:::Vector{Int}`: View of `a`.
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
