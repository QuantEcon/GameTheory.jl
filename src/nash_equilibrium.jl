"""
Finds all pure strategy Nash equilibrium for a normal form
game. It returns an empty array if there are no nash
equilibrium
"""
function pure_strategy_NE(nfg::NormalFormGame)
    # Get number of players and their actions
    np = num_players(nfg)
    na = nfg.nums_actions

    # Holder for all NE
    ne = Array(ActionProfile, 0)

    # For each action profile check whether it is NE
    as = CartesianRange(na)
    for a in as
        _a = a.I::NTuple{np, Int}
        is_nash(nfg, _a) ? push!(ne, _a) : nothing
    end

    return ne
end

