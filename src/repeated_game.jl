#=
Tools for repeated games

Authors: Chase Coleman

This file contains code to build and manage repeated games

It currently only has tools for solving two player repeated
games, but could be extended to do more with some effort.
=#

import Polyhedra: removevredundancy

"""
    RepeatedGame{N,T}

Class representing an N-player repeated game.

# Fields

- `sg::NormalFormGame{N, T}` : The stage game used to create the repeated game.
- `delta::Float64` : The common discount rate at which all players discount the
  future.
"""
struct RepeatedGame{N, T<:Real}
    sg::NormalFormGame{N, T}
    delta::Float64
end

# Type alias for 2 player game
const RepGame2 = RepeatedGame{2}

#
# Helper Functions (for 2 players)
#

"""
    RepeatedGame(p1, p2, delta)

Helper constructor that builds a repeated game for two players.

# Arguments

- `p1::Player` : The first player.
- `p2::Player` : The second player.
- `delta::Float64` : The common discount rate at which all players discount the
  future.

# Returns

- `::RepeatedGame` : The repeated game.
"""
RepeatedGame(p1::Player, p2::Player, delta::Float64) =
    RepeatedGame(NormalFormGame((p1, p2)), delta)

"""
    unpack(rpd)

Helper function that unpacks the elements of a repeated game.
"""
unpack(rpd::RepeatedGame) = (rpd.sg, rpd.delta)

# Flow utility in terms of the players actions
flow_u_1(rpd::RepGame2, a1::Int, a2::Int) =
    rpd.sg.players[1].payoff_array[a1, a2]
flow_u_2(rpd::RepGame2, a1::Int, a2::Int) =
    rpd.sg.players[2].payoff_array[a2, a1]
flow_u(rpd::RepGame2, a1::Int, a2::Int) =
    [flow_u_1(rpd, a1, a2), flow_u_2(rpd, a1, a2)]

# Computes each players best deviation given an opponent's action
best_dev_i(rpd::RepGame2, i::Int, aj::Int) =
    indmax(rpd.sg.players[i].payoff_array[:, aj])
best_dev_1(rpd::RepGame2, a2::Int) = best_dev_i(rpd, 1, a2)
best_dev_2(rpd::RepGame2, a1::Int) = best_dev_i(rpd, 2, a1)

# Computes the payoff of the best deviation
best_dev_payoff_i(rpd::RepGame2, i::Int, aj::Int) =
    maximum(rpd.sg.players[i].payoff_array[:, aj])
best_dev_payoff_1(rpd::RepGame2, a2::Int) =
    maximum(rpd.sg.players[1].payoff_array[:, a2])
best_dev_payoff_2(rpd::RepGame2, a1::Int) =
    maximum(rpd.sg.players[2].payoff_array[:, a1])

"""
    unitcircle(npts)

Places `npts` equally spaced points along the 2 dimensional circle and returns
the points with x coordinates in first column and y coordinates in second
column.
"""
function unitcircle(npts::Int)
    # Want our points placed on [0, 2π]
    incr = 2π / npts
    degrees = 0.0:incr:2π

    # Points on circle
    pts = Array{Float64}(npts, 2)
    for i=1:npts
        x = degrees[i]
        pts[i, 1] = cos(x)
        pts[i, 2] = sin(x)
    end
    return pts
end

"""
    initialize_sg_hpl(nH, o, r)

Initializes subgradients, extreme points and hyperplane levels for the
approximation of the convex value set of a 2 player repeated game.

# Arguments

- `nH::Int` : Number of subgradients used for the approximation.
- `o::Vector{Float64}` : Origin for the approximation.
- `r::Float64` : Radius for the approximation.

# Returns

- `C::Array{Float64}(nH, 1)` : The array containing the hyperplane levels.
- `H::Array{Float64}(nH, 2)` : The array containing the subgradients.
- `Z::Array{Float64}(nH, 2)` : The array containing the extreme points of the
  value set.
"""
function initialize_sg_hpl(nH::Int, o::Vector{Float64}, r::Float64)
    # First create unit circle
    H = unitcircle(nH)
    HT = H'

    # Choose origin and radius for big approximation
    Z = Array{Float64}(2, nH)
    for i=1:nH
        # We know that players can ever get worse than their
        # lowest punishment, so ignore anything below that
        Z[1, i] = o[1] + r*HT[1, i]
        Z[2, i] = o[2] + r*HT[2, i]
    end

    # Corresponding hyperplane levels
    C = squeeze(sum(HT .* Z, 1), 1)

    return C, H, Z
end

"""
    initialize_sg_hpl(rpd, nH)

Initializes subgradients, extreme points and hyperplane levels for the
approximation of the convex value set of a 2 player repeated game by choosing
an appropriate origin and radius.

# Arguments

- `rpd::RepeatedGame` : Two player repeated game.
- `nH::Int` : Number of subgradients used for the approximation.

# Returns

- `C::Array{Float64}(nH, 1)` : The array containing the hyperplane levels.
- `H::Array{Float64}(nH, 2)` : The array containing the subgradients.
- `Z::Array{Float64}(nH, 2)` : The array containing the extreme points of the
  value set.
"""
function initialize_sg_hpl(rpd::RepeatedGame, nH::Int)
    # Choose the origin to be mean of max and min payoffs
    p1_min, p1_max = extrema(rpd.sg.players[1].payoff_array)
    p2_min, p2_max = extrema(rpd.sg.players[2].payoff_array)

    o = [(p1_min + p1_max)/2.0, (p2_min + p2_max)/2.0]
    r1 = max((p1_max - o[1])^2, (o[1] - p1_min)^2)
    r2 = max((p2_max - o[2])^2, (o[2] - p2_min)^2)
    r = sqrt(r1 + r2)

    return initialize_sg_hpl(nH, o, r)
end

#
# Linear Programming Functions
#
"""
    initialize_LP_matrices(rpd, H)

Initialize matrices for the linear programming problems.

# Arguments

- `rpd::RepeatedGame` : Two player repeated game.
- `H` : The subgradients used to approximate the value set.

# Returns

- `c::Array{Float64}(nH, 1)` : Vector used to determine which subgradient is
  being used.
- `A::Array{Float64}(nH, 2)` : Matrix with nH set constraints and to be filled
  with 2 additional incentive compatibility constraints.
- `b::Array{Float64}(nH, 1)` : Vector to be filled with the values for the
  constraints.
"""
function initialize_LP_matrices(rpd::RepGame2, H)
    # Need total number of subgradients
    nH = size(H, 1)

    # Create the c vector (objective)
    c = zeros(2)

    # Create the A matrix (constraints)
    A_H = H
    A_IC_1 = zeros(1, 2)
    A_IC_2 = zeros(1, 2)
    A_IC_1[1, 1] = -rpd.delta
    A_IC_2[1, 2] = -rpd.delta
    A = vcat(A_H, A_IC_1, A_IC_2)

    # Create the b vector (constraints)
    b = Array{Float64}(nH + 2)

    return c, A, b
end

"""
    worst_value_i(rpd, H, C, i)

Given a constraint w ∈ W, this finds the worst possible payoff for agent i.

# Arugments

- `rpd::RepGame2` : Two player repeated game.
- `H::Array{Float64, 2}` : The subgradients used to approximate the value set.
- `C::Array{Float64, 1}` : The array containing the hyperplane levels.
- `i::Int` : The player of interest.
- `lp_solver` : Allows users to choose a particular solver for linear
   programming problems. Options include ClpSolver(), CbcSolver(),
   GLPKSolverLP() and GurobiSolver(). By default, it choooses ClpSolver().


# Returns

- `out::Float64` : Worst possible payoff for player i.
"""
function worst_value_i(rpd::RepGame2, H::Array{Float64, 2},
                       C::Array{Float64, 1}, i::Int,
                       lp_solver::MathProgBase.AbstractMathProgSolver=
                       ClpSolver())
    # Objective depends on which player we are minimizing
    c = zeros(2)
    c[i] = 1.0

    # Lower and upper bounds for w
    lb = [-Inf, -Inf]
    ub = [Inf, Inf]

    lpout = linprog(c, H, '<', C, lb, ub, lp_solver)
    if lpout.status == :Optimal
        out = lpout.sol[i]
    else
        out = minimum(rpd.sg.players[i].payoff_array)
    end

    return out
end

"See worst_value_i for documentation"
worst_value_1(rpd::RepGame2, H::Array{Float64, 2}, C::Array{Float64, 1},
              lp_solver::MathProgBase.AbstractMathProgSolver=ClpSolver()) =
    worst_value_i(rpd, H, C, 1, lp_solver)
"See worst_value_i for documentation"
worst_value_2(rpd::RepGame2, H::Array{Float64, 2}, C::Array{Float64, 1},
              lp_solver::MathProgBase.AbstractMathProgSolver=ClpSolver()) =
    worst_value_i(rpd, H, C, 2, lp_solver)

#
# Outer Hyper Plane Approximation
#
"""
    outerapproximation(rpd; nH=32, tol=1e-8, maxiter=500, check_pure_nash=true,
                       verbose=false, nskipprint=50,
                       plib=getlibraryfor(2, Float64),
                       lp_solver=ClpSolver())

Approximates the set of equilibrium value set for a repeated game with the
outer hyperplane approximation described by Judd, Yeltekin, Conklin 2002.

# Arguments

- `rpd::RepGame2` : Two player repeated game.
- `nH` : Number of subgradients used for the approximation.
- `tol` : Tolerance in differences of set.
- `maxiter` : Maximum number of iterations.
- `verbose` : Whether to display updates about iterations and distance.
- `nskipprint` : Number of iterations between printing information
  (assuming verbose=true).
- `check_pure_nash`: Whether to perform a check about whether a pure Nash
  equilibrium exists.
- `plib`: Allows users to choose a particular package for the geometry
  computations.
  (See [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)
  docs for more info). By default, it chooses to use SimplePolyhedraLibrary.
- `lp_solver` : Allows users to choose a particular solver for linear
   programming problems. Options include ClpSolver(), CbcSolver(),
   GLPKSolverLP() and GurobiSolver(). By default, it choooses ClpSolver().

# Returns

- `vertices::Array{Float64}` : Vertices of the outer approximation of the
  value set.
"""
function outerapproximation(rpd::RepGame2; nH=32, tol=1e-8, maxiter=500,
                            check_pure_nash=true, verbose=false, nskipprint=50,
                            plib=getlibraryfor(2, Float64),
                            lp_solver::MathProgBase.AbstractMathProgSolver=
                            ClpSolver())
    # Long unpacking of stuff
    sg, delta = unpack(rpd)
    p1, p2 = sg.players
    po_1, po_2 = p1.payoff_array, p2.payoff_array
    p1_minpayoff, p1_maxpayoff = extrema(po_1)
    p2_minpayoff, p2_maxpayoff = extrema(po_2)

    # Check to see whether at least one pure strategy NE exists
    pure_nash_exists = check_pure_nash ? length(pure_nash(sg; ntofind=1)) > 0 :
                       true
    if !pure_nash_exists
        error("No pure action Nash equilibrium exists in stage game")
    end

    # Get number of actions for each player and create action space
    nA1, nA2 = num_actions(p1), num_actions(p2)
    nAS = nA1 * nA2
    AS = QuantEcon.gridmake(1:nA1, 1:nA2)

    # Create the unit circle, points, and hyperplane levels
    C, H, Z = initialize_sg_hpl(rpd, nH)
    Cnew = copy(C)

    # Create matrices for linear programming
    c, A, b = initialize_LP_matrices(rpd, H)

    # bounds on w are [-Inf, Inf] while bounds on slack are [0, Inf]
    lb = [-Inf, -Inf]
    ub = [Inf, Inf]

    # Set iterative parameters and iterate until converged
    iter, dist = 0, 10.0
    while (iter < maxiter) & (dist > tol)
        # Compute the current worst values for each agent
        _w1 = worst_value_1(rpd, H, C, lp_solver)
        _w2 = worst_value_2(rpd, H, C, lp_solver)

        # Update all set constraints -- Copies elements 1:nH of C into b
        copy!(b, 1, C, 1, nH)

        # Iterate over all subgradients
        for ih=1:nH
            #
            # Subgradient specific instructions
            #
            h1, h2 = H[ih, :]

            # Put the right objective into c (negative because want maximize)
            c[1] = -h1
            c[2] = -h2

            # Allocate space to store all solutions
            Cia = Array{Float64}(nAS)
            Wia = Array{Float64}(2, nAS)
            for ia=1:nAS
                #
                # Action specific instruction
                #
                a1, a2 = AS[ia, :]

                # Update incentive constraints
                b[nH+1] = (1-delta)*flow_u_1(rpd, a1, a2) -
                          (1-delta)*best_dev_payoff_1(rpd, a2) - delta*_w1
                b[nH+2] = (1-delta)*flow_u_2(rpd, a1, a2) -
                          (1-delta)*best_dev_payoff_2(rpd, a1) - delta*_w2

                # Solve corresponding linear program
                lpout = linprog(c, A, '<', b, lb, ub, lp_solver)
                if lpout.status == :Optimal
                    # Pull out optimal value and compute
                    w_sol = lpout.sol
                    value = (1-delta)*flow_u(rpd, a1, a2) + delta*w_sol

                    # Save hyperplane level and continuation promises
                    Cia[ia] = h1*value[1] + h2*value[2]
                    Wia[:, ia] = value
                else
                    Cia[ia] = -Inf
                end
            end

            # Action which pushes furthest in direction h_i
            astar = indmax(Cia)
            a1star, a2star = AS[astar, :]

            # Get hyperplane level and continuation value
            Cstar = Cia[astar]
            Wstar = Wia[:, astar]
            if Cstar > -1e15
                Cnew[ih] = Cstar
            else
                error("Failed to find feasible action/continuation pair")
            end

            # Update the points
            Z[:, ih] = (1-delta)*flow_u(rpd, a1star, a2star) + delta*[Wstar[1],
                        Wstar[2]]
        end

        # Update distance and iteration counter
        dist = maximum(abs, C - Cnew)
        iter += 1

        if verbose && mod(iter, nskipprint) == 0
            println("$iter\t$dist\t($_w1, $_w2)")
        end

        if iter >= maxiter
            warn("Maximum Iteration Reached")
        end

        # Update hyperplane levels
        copy!(C, Cnew)
    end


    # Given the H-representation `(H, C)` of the computed polytope of
    # equilibrium payoff profiles, we obtain its V-representation `vertices`
    # using Polyhedra.jl (it uses `plib` which was chosen for computations)
    p = polyhedron(SimpleHRepresentation(H, C), plib)
    vertices = SimpleVRepresentation(p).V

    # Reduce the number of vertices by rounding points to the tolerance
    tol_int = round(Int, abs(log10(tol))) - 1

    # Find vertices that are unique within tolerance level
    vertices = unique(round.(vertices, tol_int), 1)

    return vertices
end


function AS{T}(rpd::RepeatedGame{2, T}; max_iter::Integer=1000,
               plib=getlibraryfor(2, T), tol::Float64=1e-5, u=nothing)

    # Initialize W0 with each entries of payoff bimatrix
    v_old = _payoff_points(rpd.sg)

    if u == nothing
        u = [minimum(rpd.sg.players[1].payoff_array),
             minimum(rpd.sg.players[2].payoff_array)]
    end

    # create VRepresentation and Polyhedron
    V = SimpleVRepresentation(v_old)
    p = polyhedron(V, plib)
    H = SimpleHRepresentation(p)
    # get rid of redundant vertices
    V = removevredundancy(V, H)
    p = polyhedron(V, plib)

    # calculate the best deviation gains
    # normalize with (1-delta)/delta
    best_dev_gains1, best_dev_gains2 = (1-rpd.delta)/rpd.delta .* _best_dev_gains(rpd.sg)

    for iter = 1:max_iter

        v_new = Vector{T}(0) # to store new vertices
        # step 1
        for a2 in 1:rpd.sg.nums_actions[2]
            for a1 in 1:rpd.sg.nums_actions[1]
                payoff1 = rpd.sg.players[1].payoff_array[a1, a2]
                payoff2 = rpd.sg.players[2].payoff_array[a2, a1]
                IC1 = u[1] + best_dev_gains1[a1, a2]
                IC2 = u[2] + best_dev_gains2[a2, a1]

                # check if the payoff point is interior
                # first check if it satisifies IC
                if all([payoff1, payoff2] .> [IC1, IC2])
                    # then check if it is in the polyhedron
                    if [payoff1, payoff2] in H
                        push!(v_new, payoff1, payoff2)
                    end
                end

                # find out the intersections of polyhedron and IC boundaries
                p_IC = polyhedron(SimpleHRepresentation(-eye(2), -[IC1, IC2]), plib)
                p_inter = intersect(p_IC, p)
                V = SimpleVRepresentation(p_inter)
                for i = 1:nvreps(V)
                    if V.V[i, 1] ≈ IC1 || V.V[i, 2] ≈ IC2
                        push!(v_new, (rpd.delta * V.V[i, :] +
                                      (1 - rpd.delta) * [payoff1, payoff2])...)
                    end
                end
            end
        end

        v_new = reshape(v_new, 2, :)'

        # get rid of redundant points
        V = SimpleVRepresentation(v_new)
        p = polyhedron(V, plib)
        H = SimpleHRepresentation(p)
        V = removevredundancy(V, H)
        p = polyhedron(V, plib)

        # check if it's converged
        # first check if the numbers of vertices are the same
        if size(v_new) == size(v_old)
            # then check the euclidean distance
            if norm(v_new-v_old) < tol
                println("converged in $(iter) iterations")
                break
            end
        end

        # check if max_iter is reached
        if iter == max_iter
            warn("Maximum Iteration Reached")
        end

        v_old = v_new

        # step 2
        # update u
        u_ = [minimum(v_new[:, 1]),
              minimum(v_new[:, 2])]
        if any(u_ .> u)
            u = u_
        end

    end

    return V.V
end

function _payoff_points{T}(g::NormalFormGame{2, T})

    nums_action_profiles = prod(g.nums_actions)
    v = Matrix{T}(nums_action_profiles, 2)
    v[:, 1] = reshape(g.players[1].payoff_array, nums_action_profiles)
    v[:, 2] = reshape(g.players[2].payoff_array', nums_action_profiles)

    return v
end

function _best_dev_gains{T}(g::NormalFormGame{2, T})

    best_dev_gains1 = (maximum(g.players[1].payoff_array, 1)
                       .- g.players[1].payoff_array)
    best_dev_gains2 = (maximum(g.players[2].payoff_array, 1)
                       .- g.players[2].payoff_array)

    return best_dev_gains1, best_dev_gains2
end
