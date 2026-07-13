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
    _check_init_pivot(init_pivot, total_num)

Check that `1 <= init_pivot <= total_num`, throwing an `ArgumentError`
otherwise.
"""
function _check_init_pivot(init_pivot::Int, total_num::Int)
    1 <= init_pivot <= total_num || throw(ArgumentError(
        "`init_pivot` must satisfy 1 <= init_pivot <= $total_num"
    ))
    return nothing
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
julia> bimatrix = [(3, 3) (3, 2)
                   (2, 2) (5, 6)
                   (0, 3) (6, 1)];

julia> g = NormalFormGame(bimatrix)
3×2 NormalFormGame{2, Int64}:
 (3, 3)  (3, 2)
 (2, 2)  (5, 6)
 (0, 3)  (6, 1)
```

Obtain a Nash equilibrium of this game by `lemke_howson` with player 1's
action 2 (out of the three actions 1, 2, and 3) as the initial pivot:

```julia
julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

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
    nums_actions = g.nums_actions
    total_num = sum(nums_actions)
    _check_init_pivot(init_pivot, total_num)
    S = float(T)

    NE = (Vector{S}(undef, nums_actions[1]), Vector{S}(undef, nums_actions[2]))
    tableaux = ntuple(i -> Matrix{S}(undef, nums_actions[3-i], total_num+1), 2)
    bases = ntuple(i -> Vector{Int}(undef, nums_actions[3-i]), 2)

    return lemke_howson!(NE, tableaux, bases, g;
                         init_pivot=init_pivot, max_iter=max_iter,
                         capping=capping, full_output=full_output)
end


"""
    lemke_howson!(NE, tableaux, bases, g; init_pivot=1, max_iter=10^6,
                  capping=nothing, full_output=Val(false),
                  col_bufs=nothing, argmins=nothing)

Same as `lemke_howson`, but allow for passing preallocated arrays `NE` (to
store the equilibrium mixed actions), `tableaux` and `bases` (for workspace),
and, as keyword arguments, the workspace arrays `col_bufs` and `argmins`.
Each keyword left as `nothing` is allocated internally, lazily.

If the players have `m` and `n` actions, `NE` must be a tuple of `Vector{S}`s
of lengths `(m, n)`, `tableaux` a tuple of `Matrix{S}`s of sizes
`(n, m+n+1)` and `(m, m+n+1)`, `bases` a tuple of `Vector{Int}`s of lengths
`(n, m)`, `col_bufs` a tuple of `Vector{S}`s of lengths `(n, m)`, and
`argmins` a `Vector{Int}` of length at least `max(m, n)`, where
`S<:AbstractFloat` (`float(T)` for the game's payoff eltype `T` under the
non-mutating `lemke_howson`).

The two members of each of `NE`, `tableaux`, and `bases` must be distinct
arrays, and `argmins` must not alias either member of `bases`; `col_bufs[1]`
and `col_bufs[2]` may be the same array when `m == n`, as each pivoting step
uses the buffer in isolation. A `DimensionMismatch` is thrown if any of the
arrays has a wrong size, and an `ArgumentError` if arrays required to be
distinct alias each other or `init_pivot` is out of range.

With `col_bufs` and `argmins` supplied, the
call performs no workspace allocations; for machine-float element types such
as `Float64` and with `full_output=Val(false)`, repeated solves then generate
no garbage-collector pressure:

    m, n = g.nums_actions
    S = Float64
    NE = (Vector{S}(undef, m), Vector{S}(undef, n))
    tableaux = (Matrix{S}(undef, n, m+n+1), Matrix{S}(undef, m, m+n+1))
    bases = (Vector{Int}(undef, n), Vector{Int}(undef, m))
    col_bufs = (Vector{S}(undef, n), Vector{S}(undef, m))
    argmins = Vector{Int}(undef, max(m, n))
    NE = lemke_howson!(NE, tableaux, bases, g;
                       col_bufs=col_bufs, argmins=argmins)
"""
function lemke_howson!(NE::NTuple{2,Vector{S}},
                       tableaux::NTuple{2,Matrix{S}},
                       bases::NTuple{2,Vector{Int}},
                       g::NormalFormGame{2,T};
                       init_pivot::Int=1,
                       max_iter::Int=10^6,
                       capping::Union{Int,Nothing}=nothing,
                       full_output::Union{Val{true},Val{false}}=Val(false),
                       col_bufs::Union{NTuple{2,Vector{S}},Nothing}=nothing,
                       argmins::Union{Vector{Int},Nothing}=nothing
                       ) where {T,S<:AbstractFloat}
    payoff_matrices = ntuple(i -> g.players[i].payoff_array, 2)
    nums_actions = g.nums_actions
    total_num = sum(nums_actions)

    _check_init_pivot(init_pivot, total_num)

    for pl in 1:2
        length(NE[pl]) == nums_actions[pl] || throw(DimensionMismatch(
            "NE[$pl] must have length $(nums_actions[pl])"))
        size(tableaux[pl]) == (nums_actions[3-pl], total_num+1) ||
            throw(DimensionMismatch(
                "tableaux[$pl] must have size ($(nums_actions[3-pl]), $(total_num+1))"))
        length(bases[pl]) == nums_actions[3-pl] || throw(DimensionMismatch(
            "bases[$pl] must have length $(nums_actions[3-pl])"))
    end
    Base.mightalias(NE[1], NE[2]) && throw(ArgumentError(
        "NE[1] and NE[2] must be separate arrays"))
    Base.mightalias(tableaux[1], tableaux[2]) && throw(ArgumentError(
        "tableaux[1] and tableaux[2] must be separate arrays"))
    Base.mightalias(bases[1], bases[2]) && throw(ArgumentError(
        "bases[1] and bases[2] must be separate arrays"))
    if col_bufs !== nothing
        for pl in 1:2
            length(col_bufs[pl]) == nums_actions[3-pl] ||
                throw(DimensionMismatch(
                    "col_bufs[$pl] must have length $(nums_actions[3-pl])"))
        end
    end
    if argmins !== nothing
        length(argmins) >= max(nums_actions...) || throw(DimensionMismatch(
            "argmins must have length at least $(max(nums_actions...))"))
        # argmins is overwritten in each pivoting step before bases[pl] is
        # read to determine the leaving variable
        for pl in 1:2
            Base.mightalias(argmins, bases[pl]) && throw(ArgumentError(
                "argmins and bases[$pl] must be separate arrays"))
        end
    end

    capping === nothing && (capping = max_iter)

    # Materialize unsupplied keyword defaults lazily
    col_bufs = col_bufs === nothing ?
        (Vector{S}(undef, nums_actions[2]), Vector{S}(undef, nums_actions[1])) :
        col_bufs
    argmins = argmins === nothing ?
        Vector{Int}(undef, max(nums_actions...)) : argmins

    converged, num_iter, init_pivot_used =
        _lemke_howson_capping!(payoff_matrices, tableaux, bases, init_pivot,
                               max_iter, capping, col_bufs, argmins)
    _get_mixed_actions!(NE, tableaux, bases)

    if full_output isa Val{false}
        return NE
    end

    res = LHResult(NE, converged, num_iter, max_iter, init_pivot_used)

    return NE, res
end



"""
    _lemke_howson_capping!(payoff_matrices, tableaux, bases, init_pivot,
                           max_iter, capping, col_bufs, argmins)

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
- `col_bufs::NTuple{2,Vector}`: Tuple of two workspace vectors of length `n`
  and `m`, respectively.
- `argmins::Vector{Int}`: Workspace vector of length at least `max(m, n)`.

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
                                capping::Int,
                                col_bufs::NTuple{2,Vector{T}},
                                argmins::Vector{Int}) where {T<:AbstractFloat}
    total = size(tableaux[2], 1) + size(tableaux[1], 1)  # m + n
    init_pivot_curr = init_pivot
    max_iter_curr = max_iter
    total_num_iter = 0

    for _ in 1:(total - 1)
        capping_curr = min(max_iter_curr, capping)

        _initialize_tableaux!(payoff_matrices, tableaux, bases)
        converged, num_iter =
            _lemke_howson_tbl!(tableaux, bases, init_pivot_curr, capping_curr,
                               col_bufs, argmins)

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
        _lemke_howson_tbl!(tableaux, bases, init_pivot_curr, max_iter_curr,
                           col_bufs, argmins)
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

    # To be added to payoffs if min <= 0
    consts = ntuple(2) do pl
        min_ = minimum(payoff_matrices[pl])
        min_ <= 0 ? convert(T, -min_ + 1) : zero(T)
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
    _lemke_howson_tbl!(tableaux, bases, init_pivot, max_iter,
                       col_bufs, argmins)

Main body of the Lemke-Howson algorithm implementation.

Perform the complementary pivoting. Modify `tableaux` and `bases` in place.

# Arguments

- `tableaux::NTuple{2,Matrix}`: Tuple of two arrays containing the tableaux,
  of shape `(n, m+n+1)` and `(m, m+n+1)`, respectively. Modified in place.
- `bases::NTuple{2,Vector{Int}}`: Tuple of two arrays containing the bases,
  of length `n` and `m`, respectively. Modified in place.
- `init_pivot::Int`: Integer `k` such that `1 <= k <= m + n`.
- `max_iter::Int`: Maximum number of pivoting steps.
- `col_bufs::NTuple{2,Vector}`: Tuple of two workspace vectors of length `n`
  and `m`, respectively.
- `argmins::Vector{Int}`: Workspace vector of length at least `max(m, n)`.

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

julia> col_bufs = (Vector{Float64}(undef, n), Vector{Float64}(undef, m));

julia> argmins = Vector{Int}(undef, max(m, n));

julia> _lemke_howson_tbl!(tableaux, bases, 2, 10, col_bufs, argmins);

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
                            max_iter::Int,
                            col_bufs::NTuple{2,Vector{T}},
                            argmins::Vector{Int}) where {T<:AbstractFloat}
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
    NE = (Vector{T}(undef, nums_actions[1]), Vector{T}(undef, nums_actions[2]))
    return _get_mixed_actions!(NE, tableaux, bases)
end

"""
    _get_mixed_actions!(NE, tableaux, bases)

In-place version of `_get_mixed_actions`: store the mixed actions in `NE`,
a tuple of vectors of lengths `m` and `n`, and return `NE`.
"""
function _get_mixed_actions!(NE::NTuple{2,Vector{T}},
                             tableaux::NTuple{2,Matrix{T}},
                             bases::NTuple{2,Vector{Int}}) where T
    nums_actions = (size(tableaux[2], 1), size(tableaux[1], 1))

    @inbounds for pl in 1:2
        out = NE[pl]
        fill!(out, zero(T))
        offset = pl == 1 ? 0 : nums_actions[1]
        start, stop = offset + 1, offset + nums_actions[pl]
        sum_ = zero(T)
        for i in 1:nums_actions[3-pl]
            k = bases[pl][i]
            if start <= k <= stop
                v = tableaux[pl][i, end]
                out[k-offset] = v
                sum_ += v
            end
        end
        if !iszero(sum_)
            out ./= sum_
        end
    end

    return NE
end
