module Games

# package code goes here
include("normal_form_game.jl")

 export Player, NormalFormGame, 
 	best_response, best_responses, is_best_response, payoff_vector,
        is_nash, pure2mixed

end # module
