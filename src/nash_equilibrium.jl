"""
Finds all pure action Nash equilibria for a normal form
game. It returns an empty array if there is no pure
action Nash.

Currently uses a brute force algorithm, but that hopefully
will change in the future.
"""
function pure_nash(nfg::NormalFormGame)
    # Get number of players and their actions
    np = num_players(nfg)
    na = nfg.nums_actions

    # Holder for all NE
    ne = Array(PureActionProfile, 0)

    # For each action profile check whether it is NE
    as = CartesianRange(na)
    for a in as
        _a = a.I
        is_nash(nfg, _a) ? push!(ne, _a) : nothing
    end

    return ne
end

