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
    "text": "Games.jl is a Julia package about algorithms and data structures for Game Theory."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Games.jl is an unregistered package that is currently under development.To install the package, open a Julia session and typePkg.clone(\"https://github.com/QuantEcon/Games.jl\")"
},

{
    "location": "index.html#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "Once installed, the Games package can be used by typingusing GamesThe Base type Player can be created by passing a payoff matrix.player1 = Player([3 1; 0 2])A 2-player NormalFormGame can be created either by passing Player instances,player2 = Player([2 0; 1 3])\ng = NormalFormGame((player1, player2))or by passing a payoff matrix directly.payoff_bimatrix = Array{Int}(2, 2, 2)\npayoff_bimatrix[1, 1, :] = [3, 2]\npayoff_bimatrix[1, 2, :] = [1, 1]\npayoff_bimatrix[2, 1, :] = [0, 0]\npayoff_bimatrix[2, 2, :] = [2, 3]\ng = NormalFormGame(payoff_bimatrix)After constructing a NormalFormGame, we can find its Nash Equilibria by using methods of Games. For example, pure_nash finds all pure action Nash Equilibria by enumeration.pure_nash(g)Please see the notebooks on QuantEcon for more details."
},

{
    "location": "index.html#notebooks-1",
    "page": "Home",
    "title": "Notebooks",
    "category": "section",
    "text": "Some notebooks for this package are available on QuantEcon. See:Tools for Game Theory\nA Recursive Formulation of Repeated Games"
},

{
    "location": "index.html#Library-Outline-1",
    "page": "Home",
    "title": "Library Outline",
    "category": "section",
    "text": "Base Types and Methods\nGame Generators\nComputing Nash Equilibria\nRepeated Games"
},

{
    "location": "lib/base_types_and_methods.html#",
    "page": "Base Types and Methods",
    "title": "Base Types and Methods",
    "category": "page",
    "text": ""
},

{
    "location": "lib/base_types_and_methods.html#base_types_and_methods-1",
    "page": "Base Types and Methods",
    "title": "Base Types and Methods",
    "category": "section",
    "text": ""
},

{
    "location": "lib/base_types_and_methods.html#Games.Action",
    "page": "Base Types and Methods",
    "title": "Games.Action",
    "category": "Constant",
    "text": "Action{T}\n\nAlias for Union{PureAction,MixedAction{T}} where T<:Real.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.ActionProfile",
    "page": "Base Types and Methods",
    "title": "Games.ActionProfile",
    "category": "Constant",
    "text": "ActionProfile\n\nAlias for Union{PureActionProfile,MixedActionProfile}.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.MixedAction",
    "page": "Base Types and Methods",
    "title": "Games.MixedAction",
    "category": "Type",
    "text": "MixedAction{T}\n\nAlias for Vector{T} where T<:Real.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.PureAction",
    "page": "Base Types and Methods",
    "title": "Games.PureAction",
    "category": "Type",
    "text": "PureAction\n\nAlias for Integer.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Type",
    "text": "NormalFormGame{N,T}\n\nType representing an N-player normal form game.\n\nFields\n\nplayers::NTuple{N,Player{N,T<:Real}} : Tuple of Player instances.\nnums_actions::NTuple{N,Int} : Tuple of the numbers of actions, one for each player.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{Array{Games.Player{N,T},1}}, Tuple{N}, Tuple{T}} where T where N",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame(players)\n\nConstructor of an N-player NormalFormGame with a vector of N Player instances.\n\nArguments\n\nplayers::Vector{Player} : Vector of Player instances.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{Array{T,2}}, Tuple{T}} where T<:Real",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame(payoffs)\n\nConstruct a symmetric 2-player NormalFormGame with a square matrix.\n\nArguments\n\npayoffs::Matrix{T<:Real} : Square matrix representing each player\'s payoff matrix.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{Array{T,M}}, Tuple{M}, Tuple{T}} where M where T<:Real",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame(payoffs)\n\nConstruct an N-player NormalFormGame for N>=2 with an array payoffs of M=N+1 dimensions, where payoffs[a_1, a_2, ..., a_N, :] contains a profile of N payoff values.\n\nArguments\n\npayoffs::Array{T<:Real} : Array with ndims=N+1 containing payoff profiles.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{Tuple{Vararg{Games.Player{N,T},N}}}, Tuple{T}} where T where N",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame(players)\n\nConstructor of an N-player NormalFormGame with a tuple of N Player instances.\n\nArguments\n\nplayers::NTuple{N,Player} : Tuple of Player instances.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{Type,Tuple{Vararg{Int64,N}}}} where N",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame([T=Float64], nums_actions)\n\nConstructor of an N-player NormalFormGame, consisting of payoffs all 0.\n\nArguments\n\nT::Type : Type of payoff values; defaults to Float64 if not specified.\nnums_actions::NTuple{N,Int} : Numbers of actions of the N players.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.NormalFormGame-Union{Tuple{N}, Tuple{T}, Tuple{Vararg{Games.Player{N,T},N} where N}} where T where N",
    "page": "Base Types and Methods",
    "title": "Games.NormalFormGame",
    "category": "Method",
    "text": "NormalFormGame(players...)\n\nConstructor of an N-player NormalFormGame with N Player instances.\n\nArguments\n\nplayers::Player{N,T}... : N Player instances\n\nExamples\n\n# p1, p2, and p3 are all of type `Player{3,T}` for some `T`\nNormalFormGame(p1, p2, p3)\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.Player",
    "page": "Base Types and Methods",
    "title": "Games.Player",
    "category": "Type",
    "text": "Player{N,T}\n\nType representing a player in an N-player normal form game.\n\nFields\n\npayoff_array::Array{T<:Real} : Array representing the player\'s payoff function.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.best_response-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void},Array{Float64,1}}",
    "page": "Base Types and Methods",
    "title": "Games.best_response",
    "category": "Method",
    "text": "best_response(player, opponents_actions, payoff_perturbation)\n\nReturn the perturbed best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1 opponents\' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent\'s mixed action) or a scalar of integer (in which case it is treated as the opponent\'s pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\npayoff_perturbation::Vector{Float64} : Vector of length equal to the number of actions of the player containing the values (\"noises\") to be added to the payoffs in determining the best response.\n\nReturns\n\n::Int : Best response action.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.best_response-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Base Types and Methods",
    "title": "Games.best_response",
    "category": "Method",
    "text": "best_response(player, opponents_actions; tie_breaking=\"smallest\", tol=1e-8)\n\nReturn a best response action to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1 opponents\' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent\'s mixed action) or a scalar of integer (in which case it is treated as the opponent\'s pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\ntie_breaking::AbstractString(\"smallest\") : Control how to break a tie (see Returns for details).\ntol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Int : If tie_breaking=\"smallest\", returns the best response action with the smallest index; if tie_breaking=\"random\", returns an action randomly chosen from the best response actions.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.best_responses-Tuple{Games.Player,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Base Types and Methods",
    "title": "Games.best_responses",
    "category": "Method",
    "text": "best_responses(player, opponents_actions; tol=1e-8)\n\nReturn all the best response actions to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1 opponents\' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent\'s mixed action) or a scalar of integer (in which case it is treated as the opponent\'s pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\nbest_responses::Vector{Int} : Vector containing all the best response actions.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_best_response-Tuple{Games.Player,Array{T,1} where T<:Real,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Base Types and Methods",
    "title": "Games.is_best_response",
    "category": "Method",
    "text": "is_best_response(player, own_action, opponents_actions; tol=1e-8)\n\nReturn true if own_action is a best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nown_action::MixedAction : Own mixed action (vector of reals).\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1 opponents\' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent\'s mixed action) or a scalar of integer (in which case it is treated as the opponent\'s pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Bool : True if own_action is a best response to opponents_actions; false otherwise.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_best_response-Tuple{Games.Player,Integer,Union{Tuple{Vararg{Array{T,1},N}} where N where T<:Real, Tuple{Vararg{T,N}} where T<:Integer where N, Union{Array{T,1}, Integer} where T<:Real, Void}}",
    "page": "Base Types and Methods",
    "title": "Games.is_best_response",
    "category": "Method",
    "text": "is_best_response(player, own_action, opponents_actions; tol=1e-8)\n\nReturn True if own_action is a best response to opponents_actions.\n\nArguments\n\nplayer::Player : Player instance.\nown_action::PureAction : Own pure action (integer).\nopponents_actions::Union{Action,ActionProfile,Void} : Profile of N-1 opponents\' actions. If N=2, then it must be a vector of reals (in which case it is treated as the opponent\'s mixed action) or a scalar of integer (in which case it is treated as the opponent\'s pure action). If N>2, then it must be a tuple of N-1 integers (pure actions) or N-1 vectors of reals (mixed actions). (For the degenerate case N=1, it must be nothing.)\n;tol::Float64 : Tolerance to be used to determine best response actions.\n\nReturns\n\n::Bool : True if own_action is a best response to opponents_actions; false otherwise.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_nash",
    "page": "Base Types and Methods",
    "title": "Games.is_nash",
    "category": "Function",
    "text": "is_nash(g, action_profile)\n\nReturn true if action_profile is a Nash equilibrium.\n\nArguments\n\ng::NormalFormGame : Instance of N-player NormalFormGame.\naction_profile::ActionProfile : Tuple of N integers (pure actions) or N vectors of reals (mixed actions).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_nash-Tuple{Games.NormalFormGame{1,T} where T<:Real,Union{Array{T,1}, Integer} where T<:Real}",
    "page": "Base Types and Methods",
    "title": "Games.is_nash",
    "category": "Method",
    "text": "is_nash(g, action)\n\nReturn true if action is a Nash equilibrium of a trivial game with 1 player.\n\nArguments\n\ng::NormalFormGame : Instance of 1-player NormalFormGame.\naction::Action : Integer (pure action) or vector of reals (mixed action).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_pareto_dominant",
    "page": "Base Types and Methods",
    "title": "Games.is_pareto_dominant",
    "category": "Function",
    "text": "is_pareto_dominant(g, action_profile)\n\nReturn true if action_profile is Pareto dominant for game g.\n\nArguments\n\ng::NormalFormGame : Instance of N-player NormalFormGame.\naction_profile::PureActionProfile : Tuple of N integers (pure actions).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.is_pareto_efficient",
    "page": "Base Types and Methods",
    "title": "Games.is_pareto_efficient",
    "category": "Function",
    "text": "is_pareto_efficient(g, action_profile)\n\nReturn true if action_profile is Pareto efficient for game g.\n\nArguments\n\ng::NormalFormGame : Instance of N-player NormalFormGame.\naction_profile::PureActionProfile : Tuple of N integers (pure actions).\n\nReturns\n\n::Bool\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.payoff_vector-Tuple{Games.Player,Tuple{Vararg{T,N}} where T<:Integer where N}",
    "page": "Base Types and Methods",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "payoff_vector(player, opponents_actions)\n\nReturn a vector of payoff values for a Player in an N>2 player game, one for each own action, given a tuple of the opponents\' pure actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::PureActionProfile : Tuple of N-1 opponents\' pure actions.\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.payoff_vector-Tuple{Games.Player{1,T} where T<:Real,Void}",
    "page": "Base Types and Methods",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "payoff_vector(player, opponent_action)\n\nReturn a vector of payoff values for a Player in a trivial game with 1 player, one for each own action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::Void\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.payoff_vector-Tuple{Games.Player{2,T} where T<:Real,Array{T,1} where T<:Real}",
    "page": "Base Types and Methods",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "payoff_vector(player, opponent_action)\n\nReturn a vector of payoff values for a Player in a 2-player game, one for each own action, given the opponent\'s mixed action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::MixedAction : Opponent\'s mixed action (vector of reals).\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.payoff_vector-Tuple{Games.Player{2,T} where T<:Real,Integer}",
    "page": "Base Types and Methods",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "payoff_vector(player, opponent_action)\n\nReturn a vector of payoff values for a Player in a 2-player game, one for each own action, given the opponent\'s pure action.\n\nArguments\n\nplayer::Player : Player instance.\nopponent_action::PureAction : Opponent\'s pure action (integer).\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.payoff_vector-Union{Tuple{Games.Player{N,T1},Tuple{Vararg{Array{T2,1},N}} where N}, Tuple{N}, Tuple{T1}, Tuple{T2}} where T2 where T1 where N",
    "page": "Base Types and Methods",
    "title": "Games.payoff_vector",
    "category": "Method",
    "text": "payoff_vector(player, opponents_actions)\n\nReturn a vector of payoff values for a Player in an N>2 player game, one for each own action, given a tuple of the opponents\' mixed actions.\n\nArguments\n\nplayer::Player : Player instance.\nopponents_actions::MixedActionProfile : Tuple of N-1 opponents\' mixed actions.\n\nReturns\n\n::Vector : Payoff vector.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.pure2mixed-Tuple{Integer,Integer}",
    "page": "Base Types and Methods",
    "title": "Games.pure2mixed",
    "category": "Method",
    "text": "pure2mixed(num_actions, action)\n\nConvert a pure action to the corresponding mixed action.\n\nArguments\n\nnum_actions::Integer : The number of the pure actions (= the length of a mixed action).\naction::PureAction : The pure action to convert to the corresponding mixed action.\n\nReturns\n\nmixed_action::Vector{Float64} : The mixed action representation of the given pure action.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Exported-1",
    "page": "Base Types and Methods",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"Games.jl\", \"normal_form_game.jl\"]\nPrivate = false"
},

{
    "location": "lib/base_types_and_methods.html#Games.MixedActionProfile",
    "page": "Base Types and Methods",
    "title": "Games.MixedActionProfile",
    "category": "Type",
    "text": "MixedActionProfile{T,N}\n\nAlias for NTuple{N,MixedAction{T}} where T<:Real.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Games.PureActionProfile",
    "page": "Base Types and Methods",
    "title": "Games.PureActionProfile",
    "category": "Type",
    "text": "PureActionProfile{N,T}\n\nAlias for NTuple{N,T} where T<:PureAction.\n\n\n\n"
},

{
    "location": "lib/base_types_and_methods.html#Internal-1",
    "page": "Base Types and Methods",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"Games.jl\", \"normal_form_game.jl\"]\nPublic = false"
},

{
    "location": "lib/game_generators.html#",
    "page": "Game Generators",
    "title": "Game Generators",
    "category": "page",
    "text": ""
},

{
    "location": "lib/game_generators.html#game_generators-1",
    "page": "Game Generators",
    "title": "Game Generators",
    "category": "section",
    "text": ""
},

{
    "location": "lib/game_generators.html#Games.covariance_game-Union{Tuple{AbstractRNG,Tuple{Vararg{Int64,N}},Real}, Tuple{N}} where N",
    "page": "Game Generators",
    "title": "Games.covariance_game",
    "category": "Method",
    "text": "covariance_game([rng=GLOBAL_RNG], nums_actions, rho)\n\nReturn a random N-player NormalFormGame instance with N>=2 where the payoff profiles are drawn independently from the standard multi-normal with the covariance of any pair of payoffs equal to rho, as studied in Rinott and Scarsini (2000).\n\nArguments\n\nrng::AbstractRNG=GLOBAL_RNG: Random number generator used.\nnums_actions::NTuple{N,Int}: Tuple of the numbers of actions, one for each player.\nrho::Real: Covariance of a pair of payoff values. Must be in [-1/(N-1), 1], where N is the number of players.\n\nReturns\n\n::NormalFormGame: The generated random N-player NormalFormGame.\n\nReferences\n\nY. Rinott and M. Scarsini, \"On the Number of Pure Strategy Nash Equilibria in Random Games,\" Games and Economic Behavior (2000), 274-293.\n\n\n\n"
},

{
    "location": "lib/game_generators.html#Games.random_game-Union{Tuple{AbstractRNG,Tuple{Vararg{Int64,N}}}, Tuple{N}} where N",
    "page": "Game Generators",
    "title": "Games.random_game",
    "category": "Method",
    "text": "random_game([rng=GLOBAL_RNG], nums_actions)\n\nReturn a random N-player NormalFormGame instance where the payoffs are drawn independently from the uniform distribution on [0, 1).\n\nArguments\n\nrng::AbstractRNG=GLOBAL_RNG: Random number generator used.\nnums_actions::NTuple{N,Int}: Tuple of the numbers of actions, one for each player.\n\nReturns\n\n::NormalFormGame: The generated random N-player NormalFormGame.\n\n\n\n"
},

{
    "location": "lib/game_generators.html#Exported-1",
    "page": "Game Generators",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"random.jl\"]\nPrivate = false"
},

{
    "location": "lib/game_generators.html#Internal-1",
    "page": "Game Generators",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"random.jl\"]\nPublic = false"
},

{
    "location": "lib/computing_nash_equilibria.html#",
    "page": "Computing Nash Equilibria",
    "title": "Computing Nash Equilibria",
    "category": "page",
    "text": ""
},

{
    "location": "lib/computing_nash_equilibria.html#computing_nash_equilibria-1",
    "page": "Computing Nash Equilibria",
    "title": "Computing Nash Equilibria",
    "category": "section",
    "text": ""
},

{
    "location": "lib/computing_nash_equilibria.html#Games.pure_nash-Tuple{Games.NormalFormGame}",
    "page": "Computing Nash Equilibria",
    "title": "Games.pure_nash",
    "category": "Method",
    "text": "pure_nash(nfg; ntofind=Inf)\n\nFinds all pure action Nash equilibria for a normal form game. It returns an empty array if there is no pure action Nash.\n\nCurrently uses a brute force algorithm, but that hopefully will change in the future.\n\nArguments\n\nnfg::NormalFormGame: Instance of N-player NormalFormGame.\nntofind::Inf: Maximal number of pure action Nash equilibria to be found; default is Inf.\n\nReturns\n\nne::Vector{NTuple{N,Int}}: Vector of pure action Nash equilibria.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Games.support_enumeration-Union{Tuple{Games.NormalFormGame{2,T}}, Tuple{T}} where T",
    "page": "Computing Nash Equilibria",
    "title": "Games.support_enumeration",
    "category": "Method",
    "text": "support_enumeration(g)\n\nCompute mixed-action Nash equilibria with equal support size for a 2-player normal form game by support enumeration. For a non-degenerate game input, these are all the Nash equilibria.\n\nThe algorithm checks all the equal-size support pairs; if the players have the same number n of actions, there are 2n choose n minus 1 such pairs. This should thus be used only for small games.\n\nArguments\n\ng::NormalFormGame{2,T}: 2-player NormalFormGame instance.\n\nReturns\n\n::Vector{NTuple{2,Vector{S}}}: Mixed-action Nash equilibria that are found, where S is Float if T is Int or Float, and Rational if T is Rational.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Games.support_enumeration_task-Tuple{Channel,Games.NormalFormGame{2,T} where T<:Real}",
    "page": "Computing Nash Equilibria",
    "title": "Games.support_enumeration_task",
    "category": "Method",
    "text": "support_enumeration_task(c, g)\n\nTask version of support_enumeration.\n\nArguments\n\nc::Channel: Channel to be binded with the support enumeration task.\ng::NormalFormGame{2}: 2-player NormalFormGame instance.\n\nReturns\n\n::Task: Runnable task for generating Nash equilibria.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Exported-1",
    "page": "Computing Nash Equilibria",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"pure_nash.jl\", \"support_enumeration.jl\"]\nPrivate = false"
},

{
    "location": "lib/computing_nash_equilibria.html#Games._indiff_mixed_action!-Union{Tuple{Array{T,2},Array{T,1},Array{T,1},Array{T,2} where T,Array{Int64,1},Array{Int64,1}}, Tuple{T}} where T<:Real",
    "page": "Computing Nash Equilibria",
    "title": "Games._indiff_mixed_action!",
    "category": "Method",
    "text": "_indiff_mixed_action!(A, b, out, payoff_matrix, own_supp, opp_supp)\n\nGiven a player\'s payoff matrix payoff_matrix, an array own_supp of this player\'s actions, and an array opp_supp of the opponent\'s actions, each of length k, compute the opponent\'s mixed action whose support equals opp_supp and for which the player is indifferent among the actions in own_supp, if any such exists. Return true if such a mixed action exists and actions in own_supp are indeed best responses to it, in which case the outcome is stored in out; false otherwise. Arrays A and b are used in intermediate steps.\n\nArguments\n\nA::Matrix{T}: Matrix of shape (k+1, k+1) used in intermediate steps, where T<:Real.\nb::Vector{T}: Vector of length k+1 used in intermediate steps, where T<:Real.\nout::Vector{T}: Vector of length k to store the nonzero values of the desired mixed action, where T<:Real.\npayoff_matrix::Matrix: The player\'s payoff matrix, of shape (m, n).\nown_supp::Vector{Int}: Vector containing the player\'s action indices, of length k.\nopp_supp::Vector{Int}: Vector containing the opponent\'s action indices, of length k.\n\nReturns\n\n::Bool: true if a desired mixed action exists and false otherwise.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Games._next_k_array!-Tuple{Array{Int64,1}}",
    "page": "Computing Nash Equilibria",
    "title": "Games._next_k_array!",
    "category": "Method",
    "text": "_next_k_array!(a)\n\nGiven an array a of k distinct nonnegative integers, return the next k-array in lexicographic ordering of the descending sequences of the elements. a is modified in place.\n\nArguments\n\na::Vector{Int}: Array of length k.\n\nReturns\n\n::Vector{Int}: Next k-array of a.\n\nExamples\n\njulia> n, k = 4, 2\n(4,2)\n\njulia> a = collect(1:k)\n2-element Array{Int64,1}:\n 1\n 2\n\njulia> while a[end] < n + 1\n           @show a\n           _next_k_array!(a)\n       end\na = [1,2]\na = [1,3]\na = [2,3]\na = [1,4]\na = [2,4]\na = [3,4]\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Games._next_k_combination-Tuple{Int64}",
    "page": "Computing Nash Equilibria",
    "title": "Games._next_k_combination",
    "category": "Method",
    "text": "_next_k_combination(x)\n\nFind the next k-combination, as described by an integer in binary representation with the k set bits, by \"Gosper\'s hack\".\n\nCopy-paste from en.wikipedia.org/wiki/Combinatorial_number_system\n\nArguments\n\nx::Int: Integer with k set bits.\n\nReturns\n\n::Int: Smallest integer > x with k set bits.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Games._support_enumeration_producer-Union{Tuple{Channel,Tuple{Array{T,2},Array{T,2}}}, Tuple{T}} where T<:Real",
    "page": "Computing Nash Equilibria",
    "title": "Games._support_enumeration_producer",
    "category": "Method",
    "text": "_support_enumeration_producer(c, payoff_matrices)\n\nMain body of support_enumeration_task.\n\nArguments\n\nc::Channel: Channel to be binded with the support enumeration task.\npayoff_matrices::NTuple{2,Matrix{T}}: Payoff matrices of player 1 and player 2, where T<:Real.\n\nPuts\n\nNTuple{2,Vector{S}}: Tuple of Nash equilibrium mixed actions, where S is Float if T is Int or Float, and Rational if T is Rational.\n\n\n\n"
},

{
    "location": "lib/computing_nash_equilibria.html#Internal-1",
    "page": "Computing Nash Equilibria",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"pure_nash.jl\", \"support_enumeration.jl\"]\nPublic = false"
},

{
    "location": "lib/repeated_games.html#",
    "page": "Repeated Games",
    "title": "Repeated Games",
    "category": "page",
    "text": ""
},

{
    "location": "lib/repeated_games.html#repeated_games-1",
    "page": "Repeated Games",
    "title": "Repeated Games",
    "category": "section",
    "text": ""
},

{
    "location": "lib/repeated_games.html#Games.RepeatedGame",
    "page": "Repeated Games",
    "title": "Games.RepeatedGame",
    "category": "Type",
    "text": "RepeatedGame{N,T}\n\nType representing an N-player repeated game.\n\nFields\n\nsg::NormalFormGame{N, T} : The stage game used to create the repeated game.\ndelta::Float64 : The common discount rate at which all players discount the future.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.RepeatedGame-Tuple{Games.Player,Games.Player,Float64}",
    "page": "Repeated Games",
    "title": "Games.RepeatedGame",
    "category": "Method",
    "text": "RepeatedGame(p1, p2, delta)\n\nHelper constructor that builds a repeated game for two players.\n\nArguments\n\np1::Player : The first player.\np2::Player : The second player.\ndelta::Float64 : The common discount rate at which all players discount the future.\n\nReturns\n\n::RepeatedGame : The repeated game.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.outerapproximation-Tuple{Games.RepeatedGame{2,T} where T<:Real}",
    "page": "Repeated Games",
    "title": "Games.outerapproximation",
    "category": "Method",
    "text": "outerapproximation(rpd; nH=32, tol=1e-8, maxiter=500, check_pure_nash=true,\n                   verbose=false, nskipprint=50,\n                   plib=getlibraryfor(2, Float64),\n                   lp_solver=ClpSolver())\n\nApproximates the set of equilibrium values for a repeated game with the outer hyperplane approximation described by Judd, Yeltekin, Conklin (2002).\n\nArguments\n\nrpd::RepGame2 : Two player repeated game.\nnH::Int : Number of subgradients used for the approximation.\ntol::Float64 : Tolerance in differences of set.\nmaxiter::Int : Maximum number of iterations.\ncheck_pure_nash: Whether to perform a check about whether a pure Nash equilibrium exists.\nverbose::Bool : Whether to display updates about iterations and distance.\nnskipprint::Int : Number of iterations between printing information (assuming verbose=true).\nplib::PolyhedraLibrary: Allows users to choose a particular package for the geometry computations. (See Polyhedra.jl docs for more info). By default, it chooses to use SimplePolyhedraLibrary.\nlp_solver::AbstractMathProgSolver : Allows users to choose a particular solver for linear programming problems. Options include ClpSolver(), CbcSolver(), GLPKSolverLP() and GurobiSolver(). By default, it choooses ClpSolver().\n\nReturns\n\nvertices::Matrix{Float64} : Vertices of the outer approximation of the value set.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.unpack-Tuple{Games.RepeatedGame}",
    "page": "Repeated Games",
    "title": "Games.unpack",
    "category": "Method",
    "text": "unpack(rpd)\n\nHelper function that unpacks the elements of a repeated game.\n\nArguments\n\nrpd::RepeatedGame : The repeated game.\n\nReturns\n\n::Tuple{NormalFormGame, Float64} : A tuple containing the stage game and the delta.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.worst_value_1",
    "page": "Repeated Games",
    "title": "Games.worst_value_1",
    "category": "Function",
    "text": "See worst_value_i for documentation\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.worst_value_2",
    "page": "Repeated Games",
    "title": "Games.worst_value_2",
    "category": "Function",
    "text": "See worst_value_i for documentation\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.worst_value_i",
    "page": "Repeated Games",
    "title": "Games.worst_value_i",
    "category": "Function",
    "text": "worst_value_i(rpd, H, C, i)\n\nGiven a constraint w ∈ W, this finds the worst possible payoff for agent i.\n\nArguments\n\nrpd::RepGame2 : Two player repeated game.\nH::Matrix{Float64} : Matrix of shape (nH, 2) containing the subgradients  here nH is the number of subgradients.\nC::Vector{Float64} : The array containing the hyperplane levels.\ni::Int : The player of interest.\nlp_solver::AbstractMathProgSolver : Allows users to choose a particular solver for linear programming problems. Options include ClpSolver(), CbcSolver(), GLPKSolverLP() and GurobiSolver(). By default, it choooses ClpSolver().\n\nReturns\n\nout::Float64 : Worst possible payoff for player i.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Exported-1",
    "page": "Repeated Games",
    "title": "Exported",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game_util.jl\", \"repeated_game.jl\"]\nPrivate = false"
},

{
    "location": "lib/repeated_games.html#Games.RepGame2",
    "page": "Repeated Games",
    "title": "Games.RepGame2",
    "category": "Type",
    "text": "RepGame2\n\nType representing a 2-player repeated game; alias for RepeatedGame{2}.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.initialize_LP_matrices-Tuple{Games.RepeatedGame{2,T} where T<:Real,Array{Float64,2}}",
    "page": "Repeated Games",
    "title": "Games.initialize_LP_matrices",
    "category": "Method",
    "text": "initialize_LP_matrices(rpd, H)\n\nInitialize matrices for the linear programming problems.\n\nArguments\n\nrpd::RepeatedGame : Two player repeated game.\nH::Matrix{Float64} : Matrix of shape (nH, 2) containing the subgradients  used to approximate the value set, where nH is the number of subgradients.\n\nReturns\n\nc::Vector{Float64} : Vector of length nH used to determine which subgradient should be used, where nH is the number of subgradients.\nA::Matrix{Float64} : Matrix of shape (nH+2, 2) with nH set constraints and to be filled with 2 additional incentive compatibility constraints.\nb::Vector{Float64} : Vector of length nH+2 to be filled with the values for the constraints.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.initialize_sg_hpl-Tuple{Games.RepeatedGame,Int64}",
    "page": "Repeated Games",
    "title": "Games.initialize_sg_hpl",
    "category": "Method",
    "text": "initialize_sg_hpl(rpd, nH)\n\nInitializes subgradients, extreme points and hyperplane levels for the approximation of the convex value set of a 2 player repeated game by choosing an appropriate origin and radius.\n\nArguments\n\nrpd::RepeatedGame : Two player repeated game.\nnH::Int : Number of subgradients used for the approximation.\n\nReturns\n\nC::Vector{Float64} : Vector of length nH containing the hyperplane levels.\nH::Matrix{Float64} : Matrix of shape (nH, 2) containing the subgradients.\nZ::Matrix{Float64} : Matrix of shape (nH, 2) containing the extreme points of the value set.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.initialize_sg_hpl-Tuple{Int64,Array{Float64,1},Float64}",
    "page": "Repeated Games",
    "title": "Games.initialize_sg_hpl",
    "category": "Method",
    "text": "initialize_sg_hpl(nH, o, r)\n\nInitializes subgradients, extreme points and hyperplane levels for the approximation of the convex value set of a 2 player repeated game.\n\nArguments\n\nnH::Int : Number of subgradients used for the approximation.\no::Vector{Float64} : Origin for the approximation.\nr::Float64 : Radius for the approximation.\n\nReturns\n\nC::Vector{Float64} : Vector of length nH containing the hyperplane levels.\nH::Matrix{Float64} : Matrix of shape (nH, 2) containing the subgradients.\nZ::Matrix{Float64} : Matrix of shape (nH, 2) containing the extreme points of the value set.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Games.unitcircle-Tuple{Int64}",
    "page": "Repeated Games",
    "title": "Games.unitcircle",
    "category": "Method",
    "text": "unitcircle(npts)\n\nPlaces npts equally spaced points along the 2 dimensional unit circle and returns the points with x coordinates in first column and y coordinates in second column.\n\nArguments\n\nnpts::Int : Number of points to be placed.\n\nReturns\n\npts::Matrix{Float64} : Matrix of shape (nH, 2) containing the coordinates of the points.\n\n\n\n"
},

{
    "location": "lib/repeated_games.html#Internal-1",
    "page": "Repeated Games",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [Games]\nPages   = [\"repeated_game_util.jl\", \"repeated_game.jl\"]\nPublic = false"
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
    "text": "Modules = [Games]"
},

]}
