using HomotopyContinuation

"""
    hc_solve(g; ntofind=Inf, options...)

Compute all isolated mixed-action Nash equilibria of an N-player normal form
game.

This function solves a system of polynomial equations arising from the
nonlinear complementarity problem representation of Nash eqiulibrium, by using
`HomotopyContinuation.jl`. For games with three or more players,
`support_enumeration`, which solves one small system per support profile, is
typically much faster.

# Arguments

- `g::NormalFormGame`: N-player NormalFormGame instance.
- `ntofind=Inf`: Maximal number of Nash equilibria to find.
- `options...`: Optional arguments to pass to `HomotopyContinuation.solve`. For
  example, the option `seed::UInt32` can set the random seed used during the
  computations; by default, the global random number generator is used
  without being reseeded. See the
  [documentation](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/solve/)
  for `HomotopyContinuation.solve` for details.

# Returns

 - `::Vector{NTuple{N,Vector{Float64}}}`: Vector of mixed-action Nash
   equilibria.

# Examples

Consider the 3-player 2-action game with 9 Nash equilibria in McKelvey and
McLennan (1996) "Computation of Equilibria in Finite Games":

```julia
julia> g = NormalFormGame((2, 2, 2));

julia> g[1, 1, 1] = [9, 8, 12];

julia> g[2, 2, 1] = [9, 8, 2];

julia> g[1, 2, 2] = [3, 4, 6];

julia> g[2, 1, 2] = [3, 4, 4];

julia> g
2×2×2 NormalFormGame{3, Float64}:
[:, :, 1] =
 (9.0, 8.0, 12.0)  (0.0, 0.0, 0.0)
 (0.0, 0.0, 0.0)   (9.0, 8.0, 2.0)

[:, :, 2] =
 (0.0, 0.0, 0.0)  (3.0, 4.0, 6.0)
 (3.0, 4.0, 4.0)  (0.0, 0.0, 0.0)

julia> Base.active_repl.options.iocontext[:compact] = true;  # Reduce digits to display

julia> NEs = hc_solve(g, show_progress=false)
9-element Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}:
 ([0.5, 0.5], [0.5, 0.5], [1.0, 0.0])
 ([0.25, 0.75], [1.0, -8.81621e-39], [0.25, 0.75])
 ([0.0, 1.0], [1.0, 0.0], [0.0, 1.0])
 ([0.0, 1.0], [0.0, 1.0], [1.0, 9.11575e-64])
 ([0.25, 0.75], [0.5, 0.5], [0.333333, 0.666667])
 ([1.0, 9.18355e-41], [1.0, -9.18355e-41], [1.0, 0.0])
 ([0.0, 1.0], [0.333333, 0.666667], [0.333333, 0.666667])
 ([0.5, 0.5], [0.333333, 0.666667], [0.25, 0.75])
 ([1.0, 0.0], [0.0, 1.0], [0.0, 1.0])

julia> all([is_nash(g, NE) for NE in NEs])
true
```
"""
function hc_solve(g::NormalFormGame{N}; ntofind=Inf, options...) where N
    ntofind <= 0 && return NTuple{N,Vector{Float64}}[]

    # Do not let HomotopyContinuation.solve reseed the global RNG by default
    options = merge((seed=nothing,), NamedTuple(options))

    f = construct_hc_system(g)

    stop_fn = isfinite(ntofind) ? r -> _is_nash(r, N) : _ -> false
    res = let n = 0
            HomotopyContinuation.solve(
                f,
                stop_early_cb = r -> stop_fn(r) && ((n += 1) >= ntofind);
                options...
            )::HomotopyContinuation.Result
    end

    NEs = [_get_action_profile(r, g.nums_actions)
           for r in res.path_results if _is_nash(r, N)]
    return NEs
end

function hc_solve(g::NormalFormGame{1}; args...)
    throw(ArgumentError("not implemented for 1-player games"))
end


function construct_hc_system(g::NormalFormGame{N}) where N
    na = g.nums_actions
    indptr = [0, cumsum(na)...]
    M = indptr[end]

    @var u[1:N] v[1:M] x[1:M]
    exs = Vector{Expression}(undef, 2M+N)

    for i in 1:N
        na_js = Base.tail((na[i:end]..., na[1:i-1]...)::NTuple{N,Int})
        for a_i in 1:na[i]
            exs[indptr[i]+a_i] = u[i] - sum(
                g.players[i].payoff_array[a_i, a_js] *
                prod([x[indptr[i+k <= N ? i+k : i+k-N]+a_js.I[k]]
                     for k in 1:N-1])
                for a_js in CartesianIndices(na_js)
            ) - v[indptr[i]+a_i]
        end
    end

    for i in 1:N
        for a_i in 1:na[i]
            exs[M+indptr[i]+a_i] = x[indptr[i]+a_i] * v[indptr[i]+a_i]
        end
    end

    for i in 1:N
        exs[2M+i] = 1 - sum(x[indptr[i]+a_i] for a_i in 1:na[i])
    end

    return System(exs)
end


function _is_nash(r::PathResult, N; nonneg_tol=1e-14)
    is_success(r) || return false
    is_real(r) || return false
    for i in N+1:length(r.solution)
        r.solution[i].re > -nonneg_tol || return false
    end
    return true
end


function _get_action_profile(r::PathResult, nums_actions::NTuple{N}) where N
    out = ntuple(i -> Array{Float64}(undef, nums_actions[i]), N)
    ind = N + sum(nums_actions)
    for i in 1:N
        for j in 1:nums_actions[i]
            ind += 1
            out[i][j] = r.solution[ind].re
        end
    end
    return out
end


# HCSolver for support_enumeration

"""
    HCSolver(; options...)

Solver for the systems of polynomial equations in `support_enumeration` for
N-player games, using `HomotopyContinuation.solve` with the polyhedral
homotopy.

# Arguments

- `options...`: Optional arguments to pass to `HomotopyContinuation.solve`.
  For example, the option `seed::UInt32` can set the random seed used during
  the computations; by default, the global random number generator is used
  without being reseeded. The defaults `compile=false`,
  `show_progress=false`, and `threading=false` are used unless overridden;
  note that compilation is much slower than interpreted evaluation for the
  small systems solved here. See the
  [documentation](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/solve/)
  for `HomotopyContinuation.solve` for details.

The systems are solved in double precision: the payoffs are converted to
`Float64` before the systems are constructed.

# Examples

```julia
julia> g = NormalFormGame((2, 2, 2));

julia> g[1, 1, 1] = [9, 8, 12];

julia> g[2, 2, 1] = [9, 8, 2];

julia> g[1, 2, 2] = [3, 4, 6];

julia> g[2, 1, 2] = [3, 4, 4];

julia> NEs = support_enumeration(g, HCSolver(seed=UInt32(1234)));

julia> length(NEs)
9
```
"""
struct HCSolver{O<:NamedTuple} <: AbstractSupportEnumerationSolver
    options::O
end

function HCSolver(; options...)
    defaults = (compile=false, show_progress=false, threading=false,
                seed=nothing)
    return HCSolver(merge(defaults, NamedTuple(options)))
end

"""
    _support_solutions(solver::HCSolver, g, supps, mixing_players)

Return the real nonsingular solutions of the indifference system on the
support profile `supps` (see `_support_equations`), computed by
`HomotopyContinuation.solve`. If the system has zero mixed volume, in which
case it has no isolated solution with all the free probabilities nonzero, an
empty vector is returned.
"""
function _support_solutions(solver::HCSolver, g::NormalFormGame{N},
                            supps, mixing_players) where N
    # HomotopyContinuation computes in double precision (and fails for
    # BigFloat coefficients in some cases); no copy if already Float64
    g = convert(NormalFormGame{N,Float64}, g)

    vars = Vector{Vector{Variable}}(undef, N)
    for i in mixing_players
        vars[i] = [Variable(:x, i, a) for a in supps[i][1:end-1]]
    end
    eqs = _support_equations(Expression, g, supps, mixing_players, vars)
    F = System(eqs,
               variables=reduce(vcat, (vars[i] for i in mixing_players)))
    res = try
        HomotopyContinuation.solve(F; solver.options...)::HomotopyContinuation.Result
    catch e
        # "Cannot compute a start system": zero mixed volume
        e isa OverflowError || rethrow()
        return Vector{Float64}[]
    end
    return real_solutions(res)
end
