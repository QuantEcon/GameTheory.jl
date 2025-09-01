"""
    pure_nash(nfg; ntofind=Inf, tol=1e-8)

Finds all pure action Nash equilibria for a normal form
game. It returns an empty array if there is no pure
action Nash.

Currently uses a brute force algorithm, but that hopefully
will change in the future.

# Arguments

- `nfg::NormalFormGame`: Instance of N-player NormalFormGame.
- `ntofind::Inf`: Maximal number of pure action Nash equilibria to be
  found; default is `prod(nfg.nums_actions)`.
- `tol::Real` : Tolerance to be used to determine best response actions.

# Returns
- `ne::Vector{NTuple{N,Int}}`: Vector of pure action Nash equilibria.

# Examples

A 4-player unanimity game example:

```julia
julia> g = NormalFormGame((2, 2, 2, 2));

julia> g[1, 1, 1, 1] = 3, 3, 3, 3;

julia> g[2, 2, 2, 2] = 4, 4, 4, 4;

julia> println(g)
2×2×2×2 NormalFormGame{4, Float64}:
[:, :, 1, 1] =
 [3.0, 3.0, 3.0, 3.0]  [0.0, 0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]

[:, :, 2, 1] =
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]

[:, :, 1, 2] =
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]

[:, :, 2, 2] =
 [0.0, 0.0, 0.0, 0.0]  [0.0, 0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0, 0.0]  [4.0, 4.0, 4.0, 4.0]

julia> pure_nash(g)
8-element Vector{NTuple{4, Int64}}:
 (1, 1, 1, 1)
 (2, 2, 1, 1)
 (2, 1, 2, 1)
 (1, 2, 2, 1)
 (2, 1, 1, 2)
 (1, 2, 1, 2)
 (1, 1, 2, 2)
 (2, 2, 2, 2)
```
"""
function pure_nash(nfg::NormalFormGame; ntofind=prod(nfg.nums_actions),
                   tol::Real=1e-8)
    # Get number of players and their actions
    np = num_players(nfg)
    na = nfg.nums_actions

    # Holder for all NE
    ne = Array{PureActionProfile{np,Int}}(undef, 0)

    # Create counter for how many to find
    nfound = 0

    for _a in CartesianIndices(na)
        if is_nash(nfg, _a.I, tol=tol)
            push!(ne, _a.I)
            nfound = nfound + 1
        end
        nfound >= ntofind && break
    end

    return ne
end
