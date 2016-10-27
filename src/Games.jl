module Games

using Clp
using MathProgBase
using QuantEcon

# package code goes here
include("normal_form_game.jl")
include("repeated_game_util.jl")
include("repeated_game.jl")

export Player, NormalFormGame,
 	   best_response, best_responses, is_best_response, payoff_vector,
       is_nash, pure2mixed, num_players, num_actions, num_opponents,
       # Repeated Games
       RepeatedGame, unpack, flow_u_1, flow_u_2, flow_u, best_dev_i,
       best_dev_1, best_dev_2, best_dev_payoff_i, best_dev_payoff_1,
       best_dev_payoff_2, worst_value_i, worst_value_1, worst_value_2,
       worst_values, outerapproximation

end # module

