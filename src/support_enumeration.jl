#=
Compute all mixed Nash equilibria of a 2-player (non-degenerate) normal
form game by support enumeration.

Julia version of QuantEcon.py/support_enumeration.py

References
----------
B. von Stengel, "Equilibrium Computation for Two-Player Games in
Strategic and Extensive Form," Chapter 3, N. Nisan, T. Roughgarden, E.
Tardos, and V. Vazirani eds., Algorithmic Game Theory, 2007.
=#

using LinearAlgebra: LAPACKException, SingularException
using QuantEcon: next_k_array!

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

# Examples

```julia
julia> bimatrix = [(3, 3) (3, 2)
                   (2, 2) (5, 6)
                   (0, 3) (6, 1)];

julia> g = NormalFormGame(bimatrix)
3×2 NormalFormGame{2, Int64}:
 (3, 3)  (3, 2)
 (2, 2)  (5, 6)
 (0, 3)  (6, 1)

julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> support_enumeration(g)
3-element Vector{Tuple{Vector{Float64}, Vector{Float64}}}:
 ([1.0, 0.0, 0.0], [1.0, 0.0])
 ([0.8, 0.2, 0.0], [0.666667, 0.333333])
 ([0.0, 0.333333, 0.666667], [0.333333, 0.666667])
```
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

- `c::Channel`: Channel to be bound to the support enumeration task.
- `g::NormalFormGame{2}`: 2-player NormalFormGame instance.

# Returns

- `::Task`: Runnable task for generating Nash equilibria.

# Examples

```julia
julia> bimatrix = [(3, 3) (3, 2)
                   (2, 2) (5, 6)
                   (0, 3) (6, 1)];

julia> g = NormalFormGame(bimatrix)
3×2 NormalFormGame{2, Int64}:
 (3, 3)  (3, 2)
 (2, 2)  (5, 6)
 (0, 3)  (6, 1)

julia> c = Channel{Tuple{Vector{Float64},Vector{Float64}}}(0);

julia> t = support_enumeration_task(c, g);

julia> bind(c, t); schedule(t);

julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> for NE in c
           display(NE)
       end
([1.0, 0.0, 0.0], [1.0, 0.0])
([0.8, 0.2, 0.0], [0.666667, 0.333333])
([0.0, 0.333333, 0.666667], [0.333333, 0.666667])
```
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

- `c::Channel`: Channel to be bound to the support enumeration task.
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
    flags_vecs = Tuple(BitVector(undef, n) for n in nums_actions)
    S = typeof(zero(T)/one(T))

    for k = 1:n_min
        supps = (collect(1:k), Vector{Int}(undef, k))
        actions = (Vector{S}(undef, k), Vector{S}(undef, k))
        A = Matrix{S}(undef, k+1, k+1)
        b = Vector{S}(undef, k+1)
        while supps[1][end] <= nums_actions[1]
            @inbounds for i in 1:k
                supps[2][i] = i
            end
            while supps[2][end] <= nums_actions[2]
                if _indiff_mixed_action!(A, b, flags_vecs[1], actions[2],
                                         payoff_matrices[1],
                                         supps[1], supps[2])
                    if _indiff_mixed_action!(A, b, flags_vecs[2], actions[1],
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

@inline function _solve!(A::Matrix, b::Vector)
    r = 0
    try
        b[:] = ldiv!(lu!(A), b)
    catch SingularException
        r = 1
    end
    return r
end

"""
    _indiff_mixed_action!(A, b, own_supp_flags, out,
                          payoff_matrix, own_supp, opp_supp)

Given a player's payoff matrix `payoff_matrix`, an array `own_supp`
of this player's actions, and an array `opp_supp` of the opponent's
actions, each of length k, compute the opponent's mixed action whose
support equals `opp_supp` and for which the player is indifferent
among the actions in `own_supp`, if any such exists. Return `true`
if such a mixed action exists and actions in `own_supp` are indeed
best responses to it, in which case the outcome is stored in `out`;
`false` otherwise. Arrays `A`, `b`, `own_supp_flags` are used in intermediate
steps.

# Arguments

- `A::Matrix{T}`: Matrix of shape (k+1, k+1) used in intermediate steps, where
  `T<:Real`.
- `b::Vector{T}`: Vector of length k+1 used in intermediate steps, where
  `T<:Real`.
- `own_supp_flags::BitVector`: BitVector of length m used in intermediate
  steps.
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
                               own_supp_flags::BitVector,
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

    own_supp_flags[:] .= false
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


# N-player support enumeration

"""
    AbstractSupportSolver

Abstract type for solvers of the systems of polynomial equations that arise in
`support_enumeration` for N-player games.

A concrete subtype `S` must implement
`_support_solutions(solver::S, g, supps, mixing)`, which returns the real
nonsingular solutions of the indifference system on the support profile
`supps` as vectors of the free probabilities (see `_support_equations`).
"""
abstract type AbstractSupportSolver end

"""
    support_enumeration(g[, solver]; ntofind=Inf, tol=1e-8)

Compute all Nash equilibria of an N-player normal form game with `N >= 3` by
support enumeration, or of a 2-player game if `solver` is given explicitly.

For each support profile, the mixed actions that make each player indifferent
among the actions in their support are the solutions of a system of
polynomial equations. The system is solved by `solver`, and the solutions are
retained as Nash equilibria if all the probabilities on the supports are
positive and no action outside the supports is a profitable deviation, as
checked by `is_nash` with tolerance `tol`.

For a regular game in the sense of Harsanyi (1973), this function returns all
the Nash equilibria; almost all games are regular. For a non-regular game, it
returns all the pure-action Nash equilibria and those mixed-action Nash
equilibria that are nonsingular solutions of their support systems, while
equilibria that are not isolated are not (fully) returned.

The number of support profiles is `prod(2^n_i - 1)`, where `n_i` is the
number of actions of player `i`, and the running time is roughly proportional
to this number (on the order of milliseconds per support profile) plus the
total number of solutions of the support systems. Support profiles in which
only one player mixes, or in which some player has more free probabilities
than the other mixing players have in total, are skipped without solving; for
2-player games this reduces to the equal-size rule. This function is
typically much faster than `hc_solve`, which solves a single large system of
polynomial equations for which the computation of the start system is
expensive, while the total number of solution paths tracked is the same.

# Arguments

- `g::NormalFormGame{N}`: N-player NormalFormGame instance.
- `solver::AbstractSupportSolver=HCSolver()`: Solver for the polynomial
  systems.
- `ntofind=Inf`: Number of Nash equilibria to find.
- `tol::Real=1e-8`: Tolerance used to check that the probabilities on the
  supports are positive and, in `is_nash`, that the mixed actions are best
  responses.

# Returns

- `::Vector{NTuple{N,Vector{Float64}}}`: Vector of mixed-action Nash
  equilibria, ordered by increasing total support size.

# Examples

Consider the 3-player 2-action game with 9 Nash equilibria in McKelvey and
McLennan (1996) "Computation of Equilibria in Finite Games":

```julia
julia> g = NormalFormGame((2, 2, 2));

julia> g[1, 1, 1] = [9, 8, 12];

julia> g[2, 2, 1] = [9, 8, 2];

julia> g[1, 2, 2] = [3, 4, 6];

julia> g[2, 1, 2] = [3, 4, 4];

julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> NEs = support_enumeration(g)
9-element Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}:
 ([1.0, 0.0], [1.0, 0.0], [1.0, 0.0])
 ([1.0, 0.0], [0.0, 1.0], [0.0, 1.0])
 ([0.0, 1.0], [1.0, 0.0], [0.0, 1.0])
 ([0.0, 1.0], [0.0, 1.0], [1.0, 0.0])
 ([0.0, 1.0], [0.333333, 0.666667], [0.333333, 0.666667])
 ([0.25, 0.75], [1.0, 0.0], [0.25, 0.75])
 ([0.5, 0.5], [0.5, 0.5], [1.0, 0.0])
 ([0.25, 0.75], [0.5, 0.5], [0.333333, 0.666667])
 ([0.5, 0.5], [0.333333, 0.666667], [0.25, 0.75])

julia> all([is_nash(g, NE) for NE in NEs])
true
```

# References

- J. C. Harsanyi, "Oddness of the Number of Equilibrium Points: A New Proof,"
  International Journal of Game Theory 2 (1973), 235-250.
- R. D. McKelvey and A. McLennan, "Computation of Equilibria in Finite Games,"
  Handbook of Computational Economics 1 (1996), 87-142.
"""
function support_enumeration(g::NormalFormGame{N}; options...) where N
    return support_enumeration(g, HCSolver(); options...)
end

function support_enumeration(g::NormalFormGame{N},
                             solver::AbstractSupportSolver;
                             ntofind=Inf, tol::Real=1e-8) where N
    N >= 2 || throw(ArgumentError("not implemented for 1-player games"))

    nums_actions = g.nums_actions
    NEs = NTuple{N,Vector{Float64}}[]

    # Support size profiles, ordered by total support size
    size_profiles =
        vec(collect(Iterators.product(ntuple(i -> 1:nums_actions[i], N)...)))
    sort!(size_profiles, by=ks -> (sum(ks), ks))

    for ks in size_profiles
        mixing = [i for i in 1:N if ks[i] > 1]
        if length(mixing) == 1
            continue  # No isolated equilibrium
        elseif length(mixing) >= 2
            # Player i's k_i - 1 indifference equations involve the other
            # mixing players' free probabilities
            num_free = sum(ks[i] - 1 for i in mixing)
            all(2 * (ks[i] - 1) <= num_free for i in mixing) || continue
        end

        supps = ntuple(i -> collect(1:ks[i]), N)
        while true
            if isempty(mixing)
                action_profile = _support_action_profile(nums_actions, supps)
                is_nash(g, action_profile, tol=tol) &&
                    push!(NEs, action_profile)
                length(NEs) >= ntofind && return NEs
            else
                for sol in _support_solutions(solver, g, supps, mixing)
                    action_profile =
                        _support_action_profile(nums_actions, supps, sol, tol)
                    action_profile === nothing && continue
                    is_nash(g, action_profile, tol=tol) || continue
                    push!(NEs, action_profile)
                    length(NEs) >= ntofind && return NEs
                end
            end
            _next_supports!(supps, nums_actions) || break
        end
    end

    return NEs
end

"""
    _next_supports!(supps, nums_actions)

Update `supps` in place to the next support profile with the same support
sizes, in the order in which the support of the last player changes fastest.
Return `false` if `supps` was the last support profile, `true` otherwise.
"""
function _next_supports!(supps::NTuple{N,Vector{Int}},
                         nums_actions::NTuple{N,Int}) where N
    for i in N:-1:1
        next_k_array!(supps[i])
        supps[i][end] <= nums_actions[i] && return true
        supps[i] .= 1:length(supps[i])
    end
    return false
end

"""
    _support_action_profile(nums_actions, supps[, sol, tol])

Return the mixed action profile with support profile `supps` whose free
probabilities are given by `sol` (see `_support_equations`), or `nothing` if
some probability on the supports is not greater than `tol`. If `sol` is
omitted, `supps` must be a pure action profile.
"""
function _support_action_profile(nums_actions::NTuple{N,Int}, supps,
                                 sol::Vector{Float64}, tol::Real) where N
    action_profile = ntuple(i -> zeros(nums_actions[i]), N)
    idx = 0
    for i in 1:N
        k = length(supps[i])
        if k == 1
            action_profile[i][supps[i][1]] = 1.
        else
            s = 0.
            for l in 1:k-1
                idx += 1
                p = sol[idx]
                p > tol || return nothing
                action_profile[i][supps[i][l]] = p
                s += p
            end
            p = 1 - s
            p > tol || return nothing
            action_profile[i][supps[i][k]] = p
        end
    end
    return action_profile
end

function _support_action_profile(nums_actions::NTuple{N,Int}, supps) where N
    action_profile = ntuple(i -> zeros(nums_actions[i]), N)
    for i in 1:N
        action_profile[i][supps[i][1]] = 1.
    end
    return action_profile
end

"""
    _support_equations(g, supps, mixing, vars)

Return the equations of the indifference system on the support profile
`supps`, as a vector of expressions in the variables `vars`.

For each player `i` in `mixing`, the players with more than one action in
their supports, `vars[i]` is a vector of `length(supps[i]) - 1` variables
representing the probabilities on the actions `supps[i][1:end-1]`, and the
probability on `supps[i][end]` is `1 - sum(vars[i])`. Players not in `mixing`
play the pure actions `supps[i][1]`. The equations are, for each `i` in
`mixing` and each `a` in `supps[i][2:end]`, the differences between the
expected payoffs of `a` and of `supps[i][1]`. The order of the variables is
that of `vars[i]` for `i` in `mixing`; the order of the equations follows
that of the variables.
"""
function _support_equations(g::NormalFormGame{N}, supps, mixing,
                            vars::Vector{<:AbstractVector}) where N
    probs_mixing = [[vars[i]; 1 - sum(vars[i])] for i in mixing]
    V = eltype(probs_mixing[1])
    probs = Vector{Vector{V}}(undef, N)  # Left undefined for pure players
    for (l, i) in enumerate(mixing)
        probs[i] = probs_mixing[l]
    end
    eqs = V[]
    for i in mixing
        opponents = ntuple(l -> mod1(i + l, N), N - 1)
        payoffs = [_support_expected_payoff(g.players[i].payoff_array, a,
                                            opponents, supps, probs, zero(V))
                   for a in supps[i]]
        for l in 2:length(supps[i])
            push!(eqs, payoffs[l] - payoffs[1])
        end
    end
    return eqs
end

"""
    _support_expected_payoff(payoff_array, a, opponents, supps, probs, z)

Return the expected payoff of action `a`, as an expression in the opponents'
probabilities `probs` restricted to their supports `supps`, where `z` is the
zero expression used to initialize the sum.
"""
function _support_expected_payoff(payoff_array::Array{T,N}, a::Int, opponents,
                                  supps, probs, z) where {T,N}
    ex = z
    for idx in CartesianIndices(ntuple(l -> length(supps[opponents[l]]), N-1))
        acts = ntuple(l -> supps[opponents[l]][idx[l]], N-1)
        coef = payoff_array[a, acts...]
        iszero(coef) && continue
        term = nothing
        for l in 1:N-1
            j = opponents[l]
            length(supps[j]) > 1 || continue  # Pure action, probability 1
            term = term === nothing ? coef * probs[j][idx[l]] :
                                      term * probs[j][idx[l]]
        end
        ex = term === nothing ? ex + coef : ex + term
    end
    return ex
end
