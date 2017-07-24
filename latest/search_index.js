var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Games.jl-1",
    "page": "Home",
    "title": "Games.jl",
    "category": "section",
    "text": "Normal Form GameComputing Nash EquilibriaPure NashRepeated GameRepeated Game Utilities\nRepeated GameRandom"
},

{
    "location": "lib/normal_form_game.html#",
    "page": "Normal Form Game",
    "title": "Normal Form Game",
    "category": "page",
    "text": ""
},

{
    "location": "lib/normal_form_game.html#normal_form_game-1",
    "page": "Normal Form Game",
    "title": "Normal Form Game",
    "category": "section",
    "text": "This is documentation for normal_form_game.jl."
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Type",
    "text": "Class representing an N-player normal form game.\n\nFields\n\nplayers::NTuple{N,Player{N,T<:Real}} : Tuple of Player instances.\nN::Int : The number of players.\nnums_actions::NTuple{N,Int} : Tuple of the numbers of actions, one for each\n\nplayer.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{Array{Games.Player{N,T},1}}, Tuple{N}, Tuple{T}} where T where N",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "Constructor of an N-player NormalFormGame.\n\nArguments\n\nplayers::Vector{Player} : Vector of Player instances.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{Array{T,2}}, Tuple{T}} where T<:Real",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame{T<:Real}(payoffs::Matrix{T})\n\nConstruct a symmetric 2-player NormalFormGame with a square matrix.\n\nArguments\n\npayoffs::Matrix{T<:Real} : Square matrix representing each player's payoff matrix.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{Array{T,M}}, Tuple{M}, Tuple{T}} where M where T<:Real",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame{T<:Real,M}(payoffs::Array{T,M})\n\nConstruct an N-player NormalFormGame for N>=2 with an array payoffs of M=N+1 dimensions, where payoffs[a_1, a_2, ..., a_N, :] contains a profile of N payoff values.\n\nArguments\n\npayoffs::Array{T<:Real} : Array with ndims=N+1 containing payoff profiles.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{Tuple{Vararg{Games.Player{N,T},N}}}, Tuple{T}} where T where N",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "Constructor of an N-player NormalFormGame.\n\nArguments\n\nplayers::NTuple{N,Player} : Tuple of Player instances.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{Type,Tuple{Vararg{Int64,N}}}} where N",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "Constructor of an N-player NormalFormGame, consisting of payoffs all 0.\n\nArguments\n\nT::Type : Type of payoff values; defaults to Float64 if not specified.\nnums_actions::NTuple{N,Int} : Numbers of actions of the N players.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{T}, Tuple{Vararg{Games.Player{N,T},N} where N}} where T where N",
    "page": "Normal Form Game",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "Constructor of an N-player NormalFormGame.\n\nArguments\n\nplayers::Player{N,T}... : N Player instances\n\nExamples\n\n# p1, p2, and p3 are all of type `Player{3,T}` for some `T`\nNormalFormGame(p1, p2, p3)\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.Player",
    "page": "Normal Form Game",
    "title": "Games.Player",
    "category": "Type",
    "text": "Type representing a player in an N-player normal form game.\n\nArguments\n\npayoff_array::Array{T<:Real} : Array representing the player's payoff\n\nfunction.\n\nFields\n\npayoff_array::Array{T<:Real} : Array representing the player's payoff\n\nfunction.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.best_response-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void},Array{Float64,1}}",
    "page": "Normal Form Game",
    "title": "Games.best_response",
    "category": "Method",
    "text": "Return the perturbed best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1\n\nopponents' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent's mixed action) or a scalar of integer (in which case it is treated as the opponent's pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n\npayoff_perturbation::Vector{Float64} : Vector of length equal to the number\n\nof actions of the player containing the values (\"noises\") to be added to the payoffs in determining the best response.\n\nReturns\n\n::Int : Best response action.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.best_response-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Normal Form Game",
    "title": "Games.best_response",
    "category": "Method",
    "text": "Return a best response action to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1\n\nopponents' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent's mixed action) or a scalar of integer (in which case it is treated as the opponent's pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n\ntie_breaking::AbstractString(\"smallest\") : Control how to break a tie (see\n\nReturns for details).\n\ntol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Int : If tie_breaking=\"smallest\", returns the best response action with\n\nthe smallest index; if tie_breaking=\"random\", returns an action randomly chosen from the best response actions.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.best_responses-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Normal Form Game",
    "title": "Games.best_responses",
    "category": "Method",
    "text": "Return all the best response actions to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1\n\nopponents' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent's mixed action) or a scalar of integer (in which case it is treated as the opponent's pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\nbest_responses::Vector{Int} : Vector containing all the best response\n\nactions.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.is_best_response-Tuple{Games.Player,Array{T,1} where T<:Real,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Normal Form Game",
    "title": "Games.is_best_response",
    "category": "Method",
    "text": "Return true if own_action is a best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nown_action::MixedAction : Own mixed action (vector of reals).\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1\n\nopponents' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent's mixed action) or a scalar of integer (in which case it is treated as the opponent's pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Bool : True if own_action is a best response to opponents_actions;\n\nfalse otherwise.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.is_best_response-Tuple{Games.Player,Integer,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Normal Form Game",
    "title": "Games.is_best_response",
    "category": "Method",
    "text": "Return True if own_action is a best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nown_action::PureAction : Own pure action (integer).\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1\n\nopponents' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent's mixed action) or a scalar of integer (in which case it is treated as the opponent's pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Bool : True if own_action is a best response to opponents_actions;\n\nvalse otherwise.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.is_nash",
    "page": "Normal Form Game",
    "title": "Games.is_nash",
    "category": "Function",
    "text": "Return true if action_profile is a Nash equilibrium.\n\nArguments\n\ng::NormalFormGame : Instance of N-player NormalFormGame.\naction_profile::ActionProfile : Tuple of N integers (pure actions) or N\n\nvectors of reals (mixed actions).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.is_nash-Tuple{Games.NormalFormGame{1,T} where T<:Real,Union{Array{T,1}, Integer} where T<:Real}",
    "page": "Normal Form Game",
    "title": "Games.is_nash",
    "category": "Method",
    "text": "Return true if action is a Nash equilibrium of a trivial game with 1 player.\n\nArguments\n\ng::NormalFormGame : Instance of 1-player NormalFormGame.\naction::Action : Integer (pure action) or vector of reals (mixed action).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.payoff_vector-Tuple{Games.Player,Tuple{Vararg{T,N}} where T<:Integer where N}",
    "page": "Normal Form Game",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "Return a vector of payoff values for a Player in an N>2 player game, one for each own action, given a tuple of the opponents' pure actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::PureActionProfile : Tuple of N-1 opponents' pure\n\nactions.\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.payoff_vector-Tuple{Games.Player{1,T} where T<:Real,Void}",
    "page": "Normal Form Game",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "Return a vector of payoff values for a Player in a trivial game with 1 player, one for each own action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::Void\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.payoff_vector-Tuple{Games.Player{2,T} where T<:Real,Array{T,1} where T<:Real}",
    "page": "Normal Form Game",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "Return a vector of payoff values for a Player in a 2-player game, one for each own action, given the opponent's mixed action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::MixedAction : Opponent's mixed action (vector of reals).\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.payoff_vector-Tuple{Games.Player{2,T} where T<:Real,Integer}",
    "page": "Normal Form Game",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "Return a vector of payoff values for a Player in a 2-player game, one for each own action, given the opponent's pure action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::PureAction : Opponent's pure action (integer).\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.payoff_vector-Union{Tuple{Games.Player{N,T1},Tuple{Vararg{Array{T2,1},N}} where N}, Tuple{N}, Tuple{T1}, Tuple{T2}} where T2 where T1 where N",
    "page": "Normal Form Game",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "Return a vector of payoff values for a Player in an N>2 player game, one for each own action, given a tuple of the opponents' mixed actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::MixedActionProfile : Tuple of N-1 opponents' mixed\n\nactions.\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Games.pure2mixed-Tuple{Integer,Integer}",
    "page": "Normal Form Game",
    "title": "Games.pure2mixed",
    "category": "Method",
    "text": "Convert a pure action to the corresponding mixed action.\n\nArguments\n\nnum_actions::Integer : The number of the pure actions (= the length of a\n\nmixed action).\n\naction::PureAction : The pure action to convert to the corresponding mixed\n\naction.\n\nReturns\n\nmixed_action::Vector{Float64} : The mixed action representation of the\n\ngiven pure action.\n\n\n\n"
},

{
    "location": "lib/normal_form_game.html#Exported-1",
    "page": "Normal Form Game",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"normal_form_game.jl\"]\nPrivate = false"
},

{
    "location": "lib/normal_form_game.html#Internal-1",
    "page": "Normal Form Game",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"normal_form_game.jl\"]\nPublic = false"
},

{
    "location": "lib/computing_nash_equilibria.html#",
    "page": "Computing Nash Equilibria",
    "title": "Computing Nash Equilibria",
    "category": "page",
    "text": ""
},

{
    "location": "lib/computing_nash_equilibria.html#Computing-Nash-Equilibria-1",
    "page": "Computing Nash Equilibria",
    "title": "Computing Nash Equilibria",
    "category": "section",
    "text": "This is documentation for Computing Nash Equilibria."
},

{
    "location": "lib/computing_nash_equilibria.html#pure_nash-1",
    "page": "Computing Nash Equilibria",
    "title": "Pure Nash",
    "category": "section",
    "text": "Documentation for pure_nash.jl."
},

{
    "location": "lib/computing_nash_equilibria.html#Games.pure_nash-Tuple{Games.NormalFormGame}",
    "page": "Computing Nash Equilibria",
    "title": "Games.pure_nash",
    "category": "Method",
    "text": "Finds all pure action Nash equilibria for a normal form game. It returns an empty array if there is no pure action Nash.\n\nCurrently uses a brute force algorithm, but that hopefully will change in the future.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Exported-1",
    "page": "Computing Nash Equilibria",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"pure_nash.jl\"]\nPrivate = false"
},

{
    "location": "lib/computing_nash_equilibria.html#Internal-1",
    "page": "Computing Nash Equilibria",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"pure_nash.jl\"]\nPublic = false"
},

{
    "location": "lib/repeated_game.html#",
    "page": "Repeated Game",
    "title": "Repeated Game",
    "category": "page",
    "text": ""
},

{
    "location": "lib/repeated_game.html#Repeated-Game-1",
    "page": "Repeated Game",
    "title": "Repeated Game",
    "category": "section",
    "text": "This is documentation for Repeated Game."
},

{
    "location": "lib/repeated_game.html#repeated_game_util-1",
    "page": "Repeated Game",
    "title": "Repeated Game Utilities",
    "category": "section",
    "text": "Documentation for repeated_game_util.jl."
},

{
    "location": "lib/repeated_game.html#Exported-1",
    "page": "Repeated Game",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game_util.jl\"]\nPrivate = false"
},

{
    "location": "lib/repeated_game.html#Games.unitcircle-Tuple{Int64}",
    "page": "Repeated Game",
    "title": "Games.unitcircle",
    "category": "Method",
    "text": "Places npts equally spaced points along the 2 dimensional circle and returns the points with x coordinates in first column and y coordinates in second column\n\ni.e. if you wanted point i, it would be pts[i, :]\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Internal-1",
    "page": "Repeated Game",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game_util.jl\"]\nPublic = false"
},

{
    "location": "lib/repeated_game.html#repeated_game-1",
    "page": "Repeated Game",
    "title": "Repeated Game",
    "category": "section",
    "text": "Documentation for repeated_game.jl."
},

{
    "location": "lib/repeated_game.html#Games.RepeatedGame",
    "page": "Repeated Game",
    "title": "Games.RepeatedGame",
    "category": "Type",
    "text": "This is a type for a specific type of repeated games\n\nIt takes a stage game that is repeated in every period and all agents discount future at rate δ\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.RepeatedGame-Tuple{Games.Player,Games.Player,Float64}",
    "page": "Repeated Game",
    "title": "Games.RepeatedGame",
    "category": "Method",
    "text": "Helper constructor that builds game from players\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.outerapproximation-Tuple{Games.RepeatedGame{2,T} where T<:Real}",
    "page": "Repeated Game",
    "title": "Games.outerapproximation",
    "category": "Method",
    "text": "Approximates the set of equilibrium value set for a repeated game with the outer hyperplane approximation described by Judd, Yeltekin, Conklin 2002\n\nNOTE: If your code fails then it might be the case that the value set is       only the value which corresponds to the pure action nash equilibrium       or you might just need more points to be precise enough for the       Polyhedra library to be able to convert to the vertice representation.\n\nThe arguments are\n\nrpd: 2 player repeated game\n\nThe keyword arguments are\n\nnH: Number of subgradients used in approximation\ntol: Tolerance in differences of set\nmaxiter: Maximum number of iterations\nverbose: Whether to display updates about iterations and distance\nnskipprint: Number of iterations between printing information (verbose=true)\ncheck_pure_nash: Whether to perform a check about whether a pure Nash equilibrium exists\nplib: Allows users to choose a particular package for the geometry computations (See Polyhedra.jl       docs for more info). By default, it chooses to use CDDLib.jl\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.unpack-Tuple{Games.RepeatedGame}",
    "page": "Repeated Game",
    "title": "Games.unpack",
    "category": "Method",
    "text": "Unpacks the elements of a repeated game\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.worst_value_1-Tuple{Games.RepeatedGame{2,T} where T<:Real,Array{Float64,2},Array{Float64,1}}",
    "page": "Repeated Game",
    "title": "Games.worst_value_1",
    "category": "Method",
    "text": "See worst_value_i for documentation\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.worst_value_2-Tuple{Games.RepeatedGame{2,T} where T<:Real,Array{Float64,2},Array{Float64,1}}",
    "page": "Repeated Game",
    "title": "Games.worst_value_2",
    "category": "Method",
    "text": "See worst_value_i for documentation\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.worst_value_i-Tuple{Games.RepeatedGame{2,T} where T<:Real,Array{Float64,2},Array{Float64,1},Int64}",
    "page": "Repeated Game",
    "title": "Games.worst_value_i",
    "category": "Method",
    "text": "Given a constraint w ∈ W, this finds the worst possible payoff for agent i\n\nThe output of this function is used to create the values associated with incentive compatibility constraints\n\nThe arguments for this function are\n\nrpd: Two player repeated game\nH: Subgradients used to approximate value set\nC: Hyperplane levels for value set approximation\ni: Which player want worst value for\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.worst_values-Tuple{Games.RepeatedGame{2,T} where T<:Real,Array{Float64,2},Array{Float64,1}}",
    "page": "Repeated Game",
    "title": "Games.worst_values",
    "category": "Method",
    "text": "See worst_value_i for documentation\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Exported-2",
    "page": "Repeated Game",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game.jl\"]\nPrivate = false"
},

{
    "location": "lib/repeated_game.html#Games.initialize_LP_matrices-Tuple{Games.RepeatedGame{2,T} where T<:Real,Any}",
    "page": "Repeated Game",
    "title": "Games.initialize_LP_matrices",
    "category": "Method",
    "text": "Initialize matrices for the linear programming problems. It sets up the A matrix since it never changes, but only allocates space for b and c since they will be filled repeatedly on the iterations\n\nWe add nH slack variables (which will be constrained to be positive) to deal with inequalities associated with Ax leq b.\n\nmin c ⋅ x     Ax < b\n\nIn this case, the c vector will be determined by which subgradient is being used, so this function only allocates space for it.\n\nThe A matrix will be filled with nH set constraints and 2 incentive compatibility constraints. The set constraints restrain the linear programming problem to pick solutions that are in the current set of continuation values while the incentive compatibility constraints ensure the agents won't deviate.\n\nThe b vector is associated with the A matrix and gives the value for constraint.\n\nThe arguments for this function are\n\nrpd: Two player repeated game\nH: The subgradients used to approximate the value set\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.initialize_sg_hpl-Tuple{Games.RepeatedGame,Int64}",
    "page": "Repeated Game",
    "title": "Games.initialize_sg_hpl",
    "category": "Method",
    "text": "This is a function that initializes the subgradients, hyperplane levels, and extreme points of the value set by choosing an appropriate origin and radius.\n\nSee initialize_sg_hpl for more documentation\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Games.initialize_sg_hpl-Tuple{Int64,Array{Float64,1},Float64}",
    "page": "Repeated Game",
    "title": "Games.initialize_sg_hpl",
    "category": "Method",
    "text": "Initializes the following things for a 2 player repeated game.\n\nsubgradients\nextreme points of the convex set for values\nhyper plane levels\n\nThese are determined in the following way\n\nSubgradients are simply chosen from the unit circle.\nThe values for the extremum of the value set are just given by choosing points along a circle with specified origin and radius\nHyperplane levels are determined by computing the hyperplane level such that the extreme points from the circle are generated\n\nThe arguments are\n\n\n\n"
},

{
    "location": "lib/repeated_game.html#Internal-2",
    "page": "Repeated Game",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game.jl\"]\nPublic = false"
},

{
    "location": "lib/random.html#",
    "page": "Random",
    "title": "Random",
    "category": "page",
    "text": ""
},

{
    "location": "lib/random.html#random-1",
    "page": "Random",
    "title": "Random",
    "category": "section",
    "text": "This is documentation for random.jl."
},

{
    "location": "lib/random.html#Games.covariance_game-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64,N}},Real}} where N",
    "page": "Random",
    "title": "Games.covariance_game",
    "category": "Method",
    "text": "covariance_game{N}(nums_actions::NTuple{N,Int}, rho::Real)\n\nReturn a random N-player NormalFormGame instance with N>=2 where the payoff profiles are drawn independently from the standard multi-normal with the covariance of any pair of payoffs equal to rho, as studied in Rinott and Scarsini (2000).\n\nArguements\n\nnums_actions::NTuple{N,Int}: Tuple of the numbers of actions,  one for each　player.\nrho::T: Covariance of a pair of payoff values. Must be in [-1/(N-1), 1], where N is the number of players.\n\nReturns\n\n::NormalFormGame: The generated random N-player NormalFormGame.\n\nReferences\n\nY. Rinott and M. Scarsini, \"On the Number of Pure Strategy Nash Equilibria in Random Games,\" Games and Economic Behavior (2000), 274-293.\n\n\n\n"
},

{
    "location": "lib/random.html#Games.random_game-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64,N}}}} where N",
    "page": "Random",
    "title": "Games.random_game",
    "category": "Method",
    "text": "random_game{N}(nums_actions::NTuple{N,Int})\n\nReturn a random N-player NormalFormGame instance where the payoffs are drawn independently from the uniform distribution on [0, 1).\n\nArguements\n\nnums_actions::NTuple{N,Int}: Tuple of the numbers of actions, one for each player.\n\nReturns\n\n::NormalFormGame: The generated random N-player NormalFormGame.\n\n\n\n"
},

{
    "location": "lib/random.html#Exported-1",
    "page": "Random",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"random.jl\"]\nPrivate = false"
},

{
    "location": "lib/random.html#Internal-1",
    "page": "Random",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"random.jl\"]\nPublic = false"
},

{
    "location": "lib/index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "lib/index.html#Index-1",
    "page": "Index",
    "title": "Index",
    "category": "section",
    "text": "modules = [Games]"
},

]}
