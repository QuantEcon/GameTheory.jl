# Katrina Evtimova, kve216

function is_supermodular(player::Player{2})
    row_diff = diff(player.payoff_array, 1)
    col_diff = diff(row_diff, 2)
    return all(col_diff .>= 0)
end
