#=
Compute mixed Nash equilibria of a 2-player normal form game by the
Lemke-Howson algorithm.
=#

using QuantEcon: _pivoting!, _lex_min_ratio_test!

"""
    LHResult

# Fields

- `NE::NTuple{2,Vector}`: Computed Nash equilibrium.
- `converged::Bool`: Whether the routine has converged.
- `num_iter::Int`: Number of iterations.
- `max_iter::Int`: Maximum number of iterations.
- `init::Int`: Initial condition used.
"""
struct LHResult{T<:Real}
    NE::NTuple{2,Vector{T}}
    converged::Bool
    num_iter::Int
    max_iter::Int
    init::Int
end


"""
    lemke_howson(g; init_pivot=1, max_iter=10^6, capping=nothing,
                 full_output=Val(false))

Find one mixed-action Nash equilibrium of a 2-player normal-form game by
the Lemke–Howson algorithm (Lemke and Howson, 1964), implemented with
"complementary pivoting" (see, e.g., von Stengel, 2007 for details).

# Arguments

- `g::NormalFormGame{2,T}`: 2-player NormalFormGame instance.
- `init_pivot::Int`: Initial pivot, an integer `k` such that `1 <= k <= m+n`,
  where integers `1, ..., m`, and `m+1, ..., m+n` correspond to the actions
  of players 1 and 2, respectively.
- `max_iter::Int`: Maximum number of pivoting steps.
- `capping::Union{Int,Nothing}`: If supplied, the routine is executed
  with the heuristic proposed by Codenotti et al. (2008); see Notes below
  for details.
- `full_output::Union{Val{true},Val{false}}`: If `Val(false)`, only the
  computed Nash equilibrium is returned. If `Val(true)`, the return value
  is `(NE, res)`, where `NE` is the Nash equilibrium and `res` is a `LHResult`
  object.

# Returns

- `NE::NTuple{2,Vector{S}}`: Tuple of computed Nash equilibrium mixed
  actions, where the type `S` is determined by `S = float(T)`.
- `res::LHResult`: Object containing information about the computation.
  Returned only when `full_output` is `Val(true)`. See `LHResult` for details.

# Examples

Consider the following game from von Stengel (2007):

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
```

Obtain a Nash equilibrium of this game by `lemke_howson` with player 1's
action 2 (out of the three actions 1, 2, and 3) as the initial pivot:

```julia
julia> NE = lemke_howson(g, init_pivot=2)
([0.0, 0.333333, 0.666667], [0.333333, 0.666667])

julia> is_nash(g, NE)
true
```

Additional information is returned if `full_output` is set `Val(true)`:

```julia
julia> NE, res = lemke_howson(g, init_pivot=2, full_output=Val(true));

julia> res.converged  # Whether the routine has converged
true

julia> res.num_iter  # Number of pivoting steps performed
4
```

# Notes

* This routine is implemented with floating-point arithmetic and thus is
  subject to numerical instability.

* If `capping` is set to a positive integer, the routine is executed with
  the heuristic proposed by Codenotti et al. (2008):

  - For `k = init_pivot, init_pivot + 1, …, init_pivot + (m+n-2)` (wrapping
    modulo `m + n` within `1:m+n`), the Lemke-Howson algorithm is executed
    with `k` as the initial pivot and `capping` as the maximum number of
    pivoting steps.

  - Otherwise, the Lemke-Howson algorithm is executed with `init_pivot +
    (m+n-1)` (wrapping modulo `m + n` within `1:m+n`) as the initial pivot,
    with a limit `max_iter` on the total number of pivoting steps.

  According to the simulation results for *uniformly random games*, for
  medium- to large-size games this heuristic outperforms the basic
  Lemke-Howson algorithm with a fixed initial pivot, where Codenotti et al.
  suggest that `capping` be set to 10.

# References

* B. Codenotti, S. De Rossi, and M. Pagan, "An Experimental Analysis of
  Lemke-Howson Algorithm," arXiv:0811.3247, 2008.
* C. E. Lemke and J. T. Howson, "Equilibrium Points of Bimatrix Games,"
  Journal of the Society for Industrial and Applied Mathematics (1964),
  413-423.
* B. von Stengel, "Equilibrium Computation for Two-Player Games in Strategic
  and Extensive Form," Chapter 3, N. Nisan, T. Roughgarden, E. Tardos, and
  V. Vazirani eds., Algorithmic Game Theory, 2007.
"""
function lemke_howson(g::NormalFormGame{2,T};
                      init_pivot::Int=1,
                      max_iter::Int=10^6,
                      capping::Union{Int,Nothing} = nothing,
                      full_output::Union{Val{true},Val{false}}=Val(false)) where T
    payoff_matrices = ntuple(i -> g.players[i].payoff_array, 2)
    nums_actions = g.nums_actions
    total_num = sum(nums_actions)

    if !(1 <= init_pivot <= total_num)
        throw(ArgumentError("`init_pivot` must satisfy 1 <= k <= $(total_num)"))
    end

    capping === nothing && (capping = max_iter)

    S = float(T)

    tableaux = ntuple(i -> Matrix{S}(undef, nums_actions[3-i], total_num+1), 2)
    bases = ntuple(i -> Vector{Int}(undef, nums_actions[3-i]), 2)

    converged, num_iter, init_pivot_used =
        _lemke_howson_capping!(payoff_matrices, tableaux, bases, init_pivot,
                               max_iter, capping)
    NE = _get_mixed_actions(tableaux, bases)

    if full_output isa Val{false}
        return NE
    end

    res = LHResult(NE, converged, num_iter, max_iter, init_pivot_used)

    return NE, res
end



"""
    _lemke_howson_capping!(payoff_matrices, tableaux, bases, init_pivot,
                           max_iter::Int, capping)

Execute the Lemke–Howson algorithm with the heuristic proposed by
Codenotti et al.

# Arguments

- `payoff_matrices::NTuple{2,Matrix}`: Tuple of two arrays representing
  payoff matrices, of shape `(m, n)` and `(n, m)`, respectively.
- `tableaux::NTuple{2,Matrix}`: Tuple of two arrays to be used to store
  the tableaux, of shape `(n, m+n+1)` and `(m, m+n+1)`, respectively.
  Modified in place.
- `bases::NTuple{2,Vector{Int}}`: Tuple of two arrays to be used to
  store the bases, of length `n` and `m`, respectively. Modified in
  place.
- `init_pivot::Int`: Integer `k` such that `1 <= k <= m + n`.
- `max_iter::Int`: Maximum number of pivoting steps.
- `capping::Int`: Value for capping. If set equal to `max_iter`, the routine
  is equivalent to the standard Lemke–Howson algorithm.

# Returns

- `converged::Bool`: Whether the pivoting terminated before `max_iter` was
  reached.
- `total_num_iter::Int`: Total number of pivoting steps performed across runs.
- `init_pivot_curr::Int`: The initial pivot used in the final run.
"""
function _lemke_howson_capping!(payoff_matrices::NTuple{2,Matrix},
                                tableaux::NTuple{2,Matrix{T}},
                                bases::NTuple{2,Vector{Int}},
                                init_pivot::Int,
                                max_iter::Int,
                                capping::Int) where {T<:AbstractFloat}
    total = size(tableaux[2], 1) + size(tableaux[1], 1)  # m + n
    init_pivot_curr = init_pivot
    max_iter_curr = max_iter
    total_num_iter = 0

    for _ in 1:(total - 1)
        capping_curr = min(max_iter_curr, capping)

        _initialize_tableaux!(payoff_matrices, tableaux, bases)
        converged, num_iter =
            _lemke_howson_tbl!(tableaux, bases, init_pivot_curr, capping_curr)

        total_num_iter += num_iter

        if converged || total_num_iter >= max_iter
            return converged, total_num_iter, init_pivot_curr
        end

        init_pivot_curr += 1
        if init_pivot_curr > total
            init_pivot_curr -= total
        end
        max_iter_curr -= num_iter
    end

    _initialize_tableaux!(payoff_matrices, tableaux, bases)
    converged, num_iter =
        _lemke_howson_tbl!(tableaux, bases, init_pivot_curr, max_iter_curr)
    total_num_iter += num_iter

    return converged, total_num_iter, init_pivot_curr
end


"""
    _initialize_tableaux!(payoff_matrices, tableaux, bases)

Given a tuple of payoff matrices, initialize the tableau and basis
arrays in place.

For each player `i`, if `minimum(payoff_matrices[i])` is non-positive,
then stored in the tableau are payoff values incremented by
`abs(minimum(payoff_matrices[i])) + 1` (to ensure the tableau does not
have a negative entry or a column identically zero).

Suppose that players 1 and 2 have `m` and `n` actions, respectively.

* `tableaux[1]` has `n` rows and `m+n+1` columns, where columns `1:m`
  and `m+1:m+n` correspond to the non-slack and slack variables,
  respectively.

* `tableaux[2]` has `m` rows and `m+n+1` columns, where columns `1:m`
  and `m+1:m+n` correspond to the slack and non-slack variables,
  respectively.

* In each `tableaux[i]`, column `m+n+1` contains the values of the basic
  variables (which are initially `1`).

* `bases[1]` and `bases[2]` contain basic variable indices, which are
  initially `m+1:m+n` and `1:m`, respectively.

# Arguments

- `payoff_matrices::NTuple{2,Matrix}`: Tuple of two arrays representing
  payoff matrices, of shape `(m, n)` and `(n, m)`, respectively.
- `tableaux::NTuple{2,Matrix}`: Tuple of two arrays to be used to store
  the tableaux, of shape `(n, m+n+1)` and `(m, m+n+1)`, respectively.
  Modified in place.
- `bases::NTuple{2,Vector{Int}}`: Tuple of two arrays to be used to
  store the bases, of length `n` and `m`, respectively. Modified in
  place.

# Returns

- `tableaux, bases`

# Examples

```julia
julia> A = [3 3; 2 5; 0 6];

julia> B = [3 2 3; 2 6 1];

julia> m, n = size(A);

julia> tableaux = (Matrix{Float64}(undef, (n, m+n+1)),
                   Matrix{Float64}(undef, (m, m+n+1)));

julia> bases = (Vector{Int}(undef, n), Vector{Int}(undef, m));

julia> tableaux, bases = _initialize_tableaux!((A, B), tableaux, bases);

julia> tableaux[1]
2×6 Matrix{Float64}:
 3.0  2.0  3.0  1.0  0.0  1.0
 2.0  6.0  1.0  0.0  1.0  1.0

julia> tableaux[2]
3×6 Matrix{Float64}:
 1.0  0.0  0.0  4.0  4.0  1.0
 0.0  1.0  0.0  3.0  6.0  1.0
 0.0  0.0  1.0  1.0  7.0  1.0

julia> bases
([4, 5], [1, 2, 3])
```
"""
function _initialize_tableaux!(payoff_matrices::NTuple{2,Matrix},
                               tableaux::NTuple{2,Matrix{T}},
                               bases::NTuple{2,Vector{Int}}) where T
    nums_actions = size(payoff_matrices[1])

    consts = zeros(T, 2)  # To be added to payoffs if min <= 0
    for pl in 1:2
        min_ = minimum(payoff_matrices[pl])
        if min_ <= 0
            consts[pl] = -min_ + 1
        end
    end

    @inbounds for (pl, (py_start, sl_start)) in enumerate(
            ((0, nums_actions[1]), (nums_actions[1], 0))
        )
        for j in 1:nums_actions[pl]
            for i in 1:nums_actions[3-pl]
                tableaux[pl][i, py_start+j] =
                    payoff_matrices[3-pl][i, j] + consts[3-pl]
            end
        end
        for j in 1:nums_actions[3-pl]
            for i in 1:nums_actions[3-pl]
                tableaux[pl][i, sl_start+j] = 0
            end
            tableaux[pl][j, sl_start+j] = 1
        end
        for i in 1:nums_actions[3-pl]
            tableaux[pl][i, end] = 1
        end

        for i in 1:nums_actions[3-pl]
            bases[pl][i] = sl_start + i
        end
    end

    return tableaux, bases
end


"""
    _lemke_howson_tbl!(tableaux, bases, init_pivot, max_iter)

Main body of the Lemke-Howson algorithm implementation.

Perform the complementary pivoting. Modify `tableaux` and `bases` in place.

# Arguments

- `tableaux::NTuple{2,Matrix}`: Tuple of two arrays containing the tableaux,
  of shape `(n, m+n+1)` and `(m, m+n+1)`, respectively. Modified in place.
- `bases::NTuple{2,Vector{Int}}`: Tuple of two arrays containing the bases,
  of length `n` and `m`, respectively. Modified in place.
- `init_pivot::Int`: Integer `k` such that `1 <= k <= m + n`.
- `max_iter::Int`: Maximum number of pivoting steps.

# Returns

- `converged::Bool`: Whether the pivoting terminated before `max_iter` was
  reached.
- `num_iter::Int`: Number of pivoting steps performed.

# Examples

```julia
julia> A = [3 3; 2 5; 0 6];

julia> B = [3 2 3; 2 6 1];

julia> m, n = size(A);

julia> tableaux = (Matrix{Float64}(undef, (n, m+n+1)),
                   Matrix{Float64}(undef, (m, m+n+1)));

julia> bases = (Vector{Int}(undef, n), Vector{Int}(undef, m));

julia> tableaux, bases = _initialize_tableaux!((A, B), tableaux, bases);

julia> _lemke_howson_tbl!(tableaux, bases, 2, 10);

julia> tableaux[1]
2×6 Matrix{Float64}:
 0.875   0.0  1.0   0.375   -0.125   0.25
 0.1875  1.0  0.0  -0.0625   0.1875  0.125

julia> tableaux[2]
3×6 Matrix{Float64}:
 1.0  -1.6         0.8  0.0  0.0  0.2
 0.0   0.466667   -0.4  1.0  0.0  0.0666667
 0.0  -0.0666667   0.2  0.0  1.0  0.133333

julia> bases
([3, 2], [1, 4, 5])
```

The outputs indicate that in the Nash equilibrium obtained, player 1's
mixed action plays actions `3` and `2` with positive weights `0.25` and
`0.125`, while player 2's mixed action plays actions `1` and `2`
(labeled as `4` and `5`) with positive weights `0.0666667` and `0.133333`.
"""
function _lemke_howson_tbl!(tableaux::NTuple{2,Matrix{T}},
                            bases::NTuple{2,Vector{Int}},
                            init_pivot::Int,
                            max_iter::Int) where {T<:AbstractFloat}
    init_player = 1
    for k in bases[1]
        if k == init_pivot
            init_player = 2
            break
        end
    end
    pls = (init_player, 3 - init_player)

    pivot = init_pivot

    m, n = (size(tableaux[2], 1), size(tableaux[1], 1))
    slack_starts = (m+1, 1)

    # Workspaces
    col_bufs = ntuple(pl -> Vector{T}(undef, size(tableaux[pl], 1)), 2)
    argmins = Vector{Int}(undef, max(m, n))

    converged = false
    num_iter  = 0

    while true
        @inbounds for pl in pls
            # Determine the leaving variable
            _, row_min = _lex_min_ratio_test!(tableaux[pl], pivot,
                                              slack_starts[pl], argmins)

            # Pivoting step: modify tableau in place
            _pivoting!(tableaux[pl], pivot, row_min, col_bufs[pl])

            # Update the basic variables and the pivot
            bases[pl][row_min], pivot = pivot, bases[pl][row_min]

            num_iter += 1

            if pivot == init_pivot
                converged = true
                break
            end
            if num_iter >= max_iter
                break
            end
        end

        if converged || num_iter >= max_iter
            break
        end
    end

    return converged, num_iter
end


"""
    _get_mixed_actions(tableaux, bases)

From `tableaux` and `bases`, extract non-slack basic variables and
return a tuple of the corresponding, normalized mixed actions.

# Arguments

- `tableaux::NTuple{2,Matrix{T}}`: Tuple of two arrays containing the tableaux,
  of shape `(n, m+n+1)` and `(m, m+n+1)`, respectively.
- `bases::NTuple{2,Vector{Int}}`: Tuple of two arrays containing the bases,
  of length `n` and `m`, respectively.

# Returns

- `::NTuple{2,Vector{T}}`: Tuple of mixed actions as given by the
  non-slack basic variables in the tableaux.
"""
function _get_mixed_actions(tableaux::NTuple{2,Matrix{T}},
                            bases::NTuple{2,Vector{Int}}) where T
    nums_actions = (size(tableaux[2], 1), size(tableaux[1], 1))
    num = nums_actions[1] + nums_actions[2]
    out = zeros(T, num)

    @inbounds for (pl, (start, stop)) in enumerate(
            ((1, nums_actions[1]), (nums_actions[1]+1, num))
        )
        sum_ = zero(T)
        for i in 1:nums_actions[3-pl]
            k = bases[pl][i]
            if start <= k <= stop
                v = tableaux[pl][i, end]
                out[k] = v
                sum_ += v
            end
        end
        if !iszero(sum_)
            @views out[start:stop] ./= sum_
        end
    end

    return out[1:nums_actions[1]], out[nums_actions[1]+1:end]
end
