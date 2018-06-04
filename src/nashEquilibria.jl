#=
Tools for computing Nash equilibria of normal form games.

Includes: Pre-processing, Lemke-Howson, Simple Search.

Authors: Arnav Sood
=#

using NonNegLeastSquares
include("Games.jl")
include("normal_form_game.jl")

"""
Render a n-player normal form game non-degenerate by adding a suitable positive constant to each payoff (i.e., by ensuring each payoff is strictly greater than 0).

##### Arguments

- `g::NormalFormGame` : A NormalFormGame object.

##### Returns

- `g::NormalFormGame` : A NormalFormGame object.
"""
function unDegenerate!(g::NormalFormGame)
    # Get the normalizing constant.
    numPlayers = num_players(g)
    minArray = Array(Real, numPlayers)
    for i in 1:numPlayers
         minArray[i] = minimum(g.players[i].payoff_array)
    end
    normalizer = abs(minimum(minArray)) + 1

    # Alter the payoff matrices.
    for i in 1:numPlayers
        g.players[i].payoff_array += normalizer
    end
    return g
end

# This is currently very slow; i.e., each call removes one dominated strategy for one player.

"""
Simplify a n-player normal form game by the iterated removal of dominated strategies.

##### Arguments

- `g::NormalFormGame` : A NormalFormGame object.

##### Returns

- `g::NormalFormGame` : A NormalFormGame object.
"""
function unDominate!(g::NormalFormGame)
    changeFlag = true
    changesMade = false
    while changeFlag == true
        changeFlag = unDominate_Auxiliary!(g::NormalFormGame)
        changesMade = changeFlag
    end
    return changesMade
end

function unDominate_Auxiliary!(g::NormalFormGame)
    # Iterate over players, flagging strictly dominated strategies, and then eliminating.
    changesMade = false
    for i in 1:num_players(g)
        player = g.players[i]
        numActions = num_actions(player)
        payoffMatrix = player.payoff_array
        dominatedActions = Array(Integer,0)
        for j in 1:numActions
            currentAction = payoffMatrix[j,:]
            for k in symdiff(numActions,j)
                 otherAction = payoffMatrix[k,:]
                 if size(otherAction.>=currentAction) == size(currentAction)
                      g.players[i].payoff_array = g.players[i].payoff_array[symdiff(numActions,j),:]
                      changesMade = true
                      return changesMade
                 end
            end
        end
    end
end

# Auxiliary function for Lemke-Howson.
function constructTableaux(g::NormalFormGame)
    # Check params.
    if num_players(g) != 2
         error("Two-player game required.")
    end

    # Bring in primitives.
    A = g.players[1].payoff_array
    B = g.players[2].payoff_array
    m = num_actions(g.players[1])
    n = num_actions(g.players[2])

    # Create tableaux.
    tableaux = Array(Matrix,2)
    tableaux[1] = cat(2, B, eye(n), ones(n, 1))
    tableaux[2] = cat(2, eye(m), A, ones(m, 1))

    return tableaux
end

# Testing implementing using another package. Need to look at syntax for passing in maxPivots.

"""
Apply the Lemke-Howson algorithm to a normal form game.

##### Arguments

- `g::NormalFormGame` : A NormalFormGame object.
- `maxPivots::Integer` : A bound on the maximum number of iterations.

##### Returns

- `nashEquilibrium::Tuple` : A Nash equilibrium for the game, given as a tuple of MixedAction objects.

"""
function applyLH(g::NormalFormGame, maxPivots::Integer=50000)
     # Pre-process.
     unDominate!(g)
     unDegenerate!(g)
     tableaux = constructTableaux(g)

     # Solve.
     nashEquilibrium = nonneg_lsq(tableaux[1],tableaux[2];alg=:pivot)
     return nashEquilibrium
end

"""
Return an ordered list of support size tuples (x, y), according to the metric introduced in this paper: http://robotics.stanford.edu/%7Eshoham/www%20papers/GEB%20computing-nash.pdf.

##### Arguments

- `g::NormalFormGame` : A 2-player NormalFormGame object.

##### Returns

- `sizeProfiles:::Array{Tuple}` : An array of integer tuples.
"""
function generateProfileList(g::NormalFormGame)
    if num_players(g) != 2
        error("Inappropriate number of players.")
    end

    m = num_actions(g.players[1])
    n = num_actions(g.players[2])
    sizeProfiles = Array{Tuple}(0)

    # Generate list.
    for i in 1:m
        for j in 1:n
            push!(sizeProfiles, (i, j))
        end
    end

    # Define the relation by overloading the < relation on tuples.
    function isless(x::Tuple, y::Tuple)
        if abs(x[1]-x[2]) < abs(y[1]-y[2])
            return true
        elseif abs(x[1]-x[2]) == abs(y[1]-y[2]) && (x[1]+x[2]) < (y[1]+y[2])
            return true
        else
            return false
        end
    end

    sort!(sizeProfiles)
    return sizeProfiles
end

"""
Determine whether a given action for a given player is conditionally dominated, i.e., is strictly dominated given a particular restriction of the opponent's strategies. The tuple is a (player, support) tuple, where player is an integer and support is an array. The other integer indexes the protagonist's action.

##### Arguments

- `g::NormalFormGame` : A NormalFormGame object.
- `oppTuple::Tuple` : A (i, [x, y, z..]) tuple, where i indexes a player and the array lists her supported actions.
- `protAction::Integer` : The index for an action for the player not indexed by i.
"""
function isConditionallyDominated(g::NormalFormGame, oppTuple::Tuple,     protAction::Integer)

# Extract primitives.
opponent = aTuple[1]
opponentActions = aTuple[2]
protagonist = symdiff([1,2],opponent)[1]
isDominated = false

# Copy the game, and eliminate actions we've restricted away (i.e., eliminate opponent's rows, and our columns.)
g_temp = copy(g)
g_temp.players[opponent].payoff_array = g.players[opponent].payoff_array[opponentActions,:]
g_temp.players[protagonist].payoff_array = g.players[protagonist].payoff_array[:,opponentActions]

# Now, check to see whether the given action is conditionally dominated.
newPayoffs = g_temp.players[protagonist].payoff_array
numActions = num_actions(g.players[protagonist])
for 1 in 1:numActions
    if size(newPayoffs[i,:].>newPayoffs[protActions,:]) == size(newPayoffs[i,:])
        isDominated = true
        return isDominated
    end
end

return isDominated
end
