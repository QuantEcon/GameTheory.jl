module Games

# Packages

# package code goes here
include("normal_form_game.jl")
include("nash_equilibrium.jl")

export Player, NormalFormGame,  # Types

       # Type aliases
       Action, MixedAction, PureAction, ActionProfile,

       # Normal form game functions
 	   best_response, best_responses, is_best_response, payoff_vector,
       is_nash, pure2mixed, pure_strategy_NE,

       # General functions
       num_players, num_actions, num_opponents

end # module
