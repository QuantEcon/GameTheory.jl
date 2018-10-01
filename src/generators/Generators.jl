module Generators

using Random

import Games: Player, NormalFormGame

include("bimatrix_generators.jl")
export blotto_game, ranking_game, sgc_game, unit_vector_game, tournament_game

end # module
