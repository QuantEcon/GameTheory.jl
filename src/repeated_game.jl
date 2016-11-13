#=
This file contains code to build and manage repeated games

It currently only has tools for solving two player repeated
games, but could be extended to do more with some effort.
=#
using Polyhedra
using CDDLib
const _polyhedra_lib = CDDLibrary


"""
This is a type for a specific type of repeated games

It takes a stage game that is repeated in every period
and all agents discount future at rate δ
"""
immutable RepeatedGame{N, T<:Real}
    sg::NormalFormGame{N, T}
    δ::Float64
end

#
# Helper Functions (for 2 players)
#

"""
Helper constructor that builds game from players
"""
RepeatedGame(p1::Player, p2::Player, δ::Float64) =
    RepeatedGame(NormalFormGame((p1, p2)), δ)

"""
Unpacks the elements of a repeated game
"""
unpack(rpd::RepeatedGame) = (rpd.sg, rpd.δ)

# Flow utility in terms of the players actions
flow_u_1(rpd, a1::Int, a2::Int) = (1-rpd.δ)*rpd.sg.players[1].payoff_array[a1, a2]
flow_u_2(rpd, a1::Int, a2::Int) = (1-rpd.δ)*rpd.sg.players[2].payoff_array[a2, a1]
flow_u(rpd, a1::Int, a2::Int) = [flow_u_1(rpd, a1, a2), flow_u_2(rpd, a1, a2)]

# Computes each players best deviation given an opponent's action
best_dev_i(rpd, i::Int, aj::Int) = indmax(rpd.sg.players[i].payoff_array[:, aj])
best_dev_1(rpd, a2::Int) = best_dev_i(rpd, 1, a2)
best_dev_2(rpd, a1::Int) = best_dev_i(rpd, 2, a1)

# Computes the payoff of the best deviation
best_dev_payoff_i(rpd, i::Int, aj::Int) = (1-rpd.δ)*maximum(rpd.sg.players[i].payoff_array[:, aj])
best_dev_payoff_1(rpd, a2::Int) = (1-rpd.δ)*maximum(rpd.sg.players[1].payoff_array[:, a2])
best_dev_payoff_2(rpd, a1::Int) = (1-rpd.δ)*maximum(rpd.sg.players[2].payoff_array[:, a1])

"""
Creates the unit circle of continuation value points starting at origin `o`
with radius `r`. Additionally, this function crops any points that go
beyond the bounds permitted by the game -- For example, if the max payoff
for player 1 is 5 then any point that awards player 1 higher than 5
continuation utility is cropped to 5.
"""
function create_unitcircle_points{T<:Real}(rpd::RepeatedGame{2, T}, nH::Int,
                                           o::Vector{Float64}, r::Float64)
    # First create unit circle
    H = unitcircle(nH)
    HT = H'

    # Choose origin and radius for big approximation
    p1_min, p1_max = extrema(rpd.sg.players[1].payoff_array)
    p2_min, p2_max = extrema(rpd.sg.players[2].payoff_array)
    Z = Array(Float64, 2, nH)
    for i=1:nH
        # We know that players can ever get worse than their
        # lowest punishment, so ignore anything below that
        Z[1, i] = min(max(o[1] + r*HT[1, i], p1_min), p1_max)
        Z[2, i] = min(max(o[2] + r*HT[2, i], p2_min), p2_max)
    end

    # Corresponding hyperplane levels
    C = squeeze(sum(HT .* Z, 1), 1)

    return C, H, Z
end

function create_unitcircle_points(rpd::RepeatedGame, nH::Int)

    # Choose the origin to be mean of max and min payoffs
    p1_min, p1_max = extrema(rpd.sg.players[1].payoff_array)
    p2_min, p2_max = extrema(rpd.sg.players[2].payoff_array)

    o = [(p1_min + p1_max)/2.0, (p2_min + p2_max)/2.0]
    r1 = max((p1_max - o[1])^2, (o[1] - p1_min)^2)
    r2 = max((p2_max - o[2])^2, (o[2] - p2_min)^2)
    r = sqrt(r1 + r2)

    return create_unitcircle_points(rpd, nH, o, r)
end

#
# Linear Programming Functions
#
"""
Initialize matrices for the linear programming problems. It sets up the A
matrix since it never changes, but only allocates space for b and c since
they will be filled repeatedly on the iterations

We add nH slack variables (which will be constrained to be positive) to
deal with inequalities associated with Ax \leq b.

min c ⋅ x
    Ax = b

In this case, the `c` vector will be determined by which subgradient is being
used, so this function only allocates space for it.

The `A` matrix will be filled with nH set constraints and 2 incentive compatibility
constraints. The set constraints restrain the linear programming problem to pick
solutions that are in the current set of continuation values while the incentive
compatibility constraints ensure the agents won't deviate.

The `b` vector is associated with the `A` matrix and gives the value for constraint.
"""
function initialize_LP_matrices(rpd, H)
    # Need slack variable for every subgradient and additional 2 incentive constraint
    nH = size(H, 1)
    nslack = nH + 2

    # Create the c vector (objective)
    c = zeros(nslack+2)

    # Create the A matrix (constraints)
    A_H = hcat(H, eye(nH, nslack))
    A_IC_1 = zeros(1, nslack+2)
    A_IC_2 = zeros(1, nslack+2)
    A_IC_1[1] = -rpd.δ
    A_IC_1[end-1] = 1.0
    A_IC_2[2] = -rpd.δ
    A_IC_2[end] = 1.0
    A = vcat(A_H, A_IC_1, A_IC_2)

    # Create the b vector (constraints)
    b = Array(Float64, nslack)

    return c, A, b
end

"""
Given a constraint w ∈ W, this finds the worst possible payoff for agent i

The output of this function is used to create the values associated with
incentive compatibility constraints
"""
function worst_value_i(rpd::RepeatedGame, H::Array{Float64, 2}, C::Array{Float64, 1}, i::Int)
    # Objective depends on which player we are minimizing
    c = zeros(2)
    c[i] = 1.0

    # Lower and upper bounds for w
    lb = [-Inf, -Inf]
    ub = [Inf, Inf]

    lpout = linprog(c, H, '<', C, lb, ub, ClpSolver())
    if lpout.status == :Optimal
        out = lpout.sol[i]
    else
        out = minimum(rpd.sg.players[i].payoff_array)
    end

    return out
end

worst_value_1(rpd::RepeatedGame, H::Array{Float64, 2}, C::Array{Float64, 1}) =
    worst_value_i(rpd, H, C, 1)
worst_value_2(rpd::RepeatedGame, H::Array{Float64, 2}, C::Array{Float64, 1}) =
    worst_value_i(rpd, H, C, 2)
worst_values(rpd::RepeatedGame, H::Array{Float64, 2}, C::Array{Float64, 1}) =
    (worst_value_1(rpd, H, C), worst_value_2(rpd, H, C))

#
# Outer Hyper Plane Approximation
#
"""
Approximates the set of equilibrium value set for a repeated game with the
outer hyperplane approximation described by Judd, Yeltekin, Conklin 2002
"""
function outerapproximation(rpd::RepeatedGame; nH=32, tol=1e-8, maxiter=500, nskipprint=1)
    # Long unpacking of stuff
    sg, δ = unpack(rpd)
    p1, p2 = sg.players
    po_1, po_2 = p1.payoff_array, p2.payoff_array
    p1_minpayoff, p1_maxpayoff = extrema(po_1)
    p2_minpayoff, p2_maxpayoff = extrema(po_2)

    # Get number of actions for each player
    nA1, nA2 = num_actions(p1), num_actions(p2)
    nAS = nA1 * nA2

    # Create action space
    AS = QuantEcon.gridmake(1:nA1, 1:nA2)

    # Create the unit circle, points, and hyperplane levels
    C, H, Z = create_unitcircle_points(rpd, nH)
    Cnew = copy(C)

    # Create matrices for linear programming
    c, A, b = initialize_LP_matrices(rpd, H)

    # bounds on w are [-Inf, Inf] while bounds on slack are [0, Inf]
    lb = vcat([-Inf, -Inf], zeros(nH+2))
    ub = fill(Inf, nH+4)

    # Set iterative parameters and iterate until converged
    iter, dist = 0, 10.0
    while (iter < maxiter) && (dist > tol)
        # Compute the current worst values for each agent
        _w1 = worst_value_1(rpd, H, C)
        _w2 = worst_value_2(rpd, H, C)

        # Update all set constraints
        copy!(b, 1, C)

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
            Cia = Array(Float64, nAS)
            Wia = Array(Float64, 2, nAS)
            for ia=1:nAS
                #
                # Action specific instruction
                #
                a1, a2 = AS[ia, :]

                # Update incentive constraints
                b[nH+1] = flow_u_1(rpd, a1, a2) - best_dev_payoff_1(rpd, a2) - δ*_w1
                b[nH+2] = flow_u_2(rpd, a1, a2) - best_dev_payoff_2(rpd, a1) - δ*_w2

                # Solve corresponding linear program
                lpout = linprog(c, A, '=', b, lb, ub, ClpSolver())
                if lpout.status == :Optimal
                    # Pull out optimal value and compute
                    w_sol = lpout.sol[1:2]
                    value = flow_u(rpd, a1, a2) + δ*w_sol

                    # Save hyperplane level and continuation promises
                    Cia[ia] = h1*value[1] + h2*value[2]
                    Wia[:, ia] = value
                else
                    @show lpout.status
                    Cia[ia] = -Inf
                end
            end

            # Action which pushes furthest in direction h_i
            astar = indmax(Cia)
            a1star, a2star = AS[astar, :]

            # Get hyperplane level and continuation value
            Cstar = Cia[astar]
            Wstar = Wia[:, astar]
            if Cstar > -Inf
                Cnew[ih] = Cstar
            else
                println("FAIL")
                Cnew[ih] = h1*p1_minpayoff + h2*p2_minpayoff
            end

            # Update the points
            Z[:, ih] = flow_u(rpd, a1star, a2star) + δ*[Wstar[1], Wstar[2]]
        end

        # Update distance and iteration counter
        dist = maxabs(C - Cnew)
        iter += 1
        mod(iter, nskipprint) == 0 ? println("$iter\t$dist\t($_w1, $_w2)") : nothing

        # Update hyperplane levels
        copy!(C, Cnew)
    end

    # Given the H-representation `(H, C)` of the computed polytope of
    # equilibrium payoff profiles, we obtain its V-representation `vertices`.
    # Here we use CDDLib.jl, a Julia wrapper of cdd, through Polyhedra.jl.
    hrep = SimpleHRepresentation(H, C)
    poly = polyhedron(hrep, _polyhedra_lib())
    vrep = getvrep(poly)
    vertices = SimpleVRepresentation(vrep).V

    return vertices
end
