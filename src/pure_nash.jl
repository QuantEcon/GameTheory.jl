"""
    pure_nash(nfg; ntofind=Inf)

Finds all pure action Nash equilibria for a normal form
game. It returns an empty array if there is no pure
action Nash.

Currently uses a brute force algorithm, but that hopefully
will change in the future.

# Arguments

- `nfg::NormalFormGame`: Instance of N-player NormalFormGame.
- `ntofind::Inf`: Maximal number of pure nash action Nash
  Equilibria to be found; default is `Inf`.

# Returns
- `ne::Vector{NTuple{na,Int}}`: Vector of pure action Nash
  Equilibria. `na` is the number of players.
"""
function pure_nash(nfg::NormalFormGame; ntofind=Inf)
    # Get number of players and their actions
    np = num_players(nfg)
    na = nfg.nums_actions

    # Holder for all NE
    ne = Array{PureActionProfile{np,Int}}(0)

    # For each action profile check whether it is NE
    as = CartesianRange(na)

    # Create counter for how many to find and iterator for actions
    nfound = 0
    state = start(as)

    while (nfound < ntofind) & !done(as, state)
        # Get next action pair
        a, state = next(as, state)
        _a = a.I

        if is_nash(nfg, _a)
            push!(ne, _a)
            nfound = nfound + 1
        end
    end

    return ne
end
