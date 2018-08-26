module Generators

# 0.6/0.7 compatibility
using Compat
using Compat.Random
import Compat.cumsum!
import Compat.maximum
import Compat.sum

import Games: Player, NormalFormGame

include("bimatrix_generators.jl")
export blotto_game, ranking_game, sgc_game, unit_vector_game, tournament_game

end # module
