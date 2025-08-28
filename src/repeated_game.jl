#=
Tools for repeated games

This file contains code to build and manage repeated games

It currently only has tools for solving two player repeated
games, but could be extended to do more with some effort.
=#

using Polyhedra
using MathOptInterface
using HiGHS

const MOI = MathOptInterface

"""
    RepeatedGame{N,T,TD}

Type representing an N-player repeated game.

# Fields

- `sg::NormalFormGame{N, T}` : The stage game used to create the repeated game.
- `delta::TD` : The common discount rate at which all players discount the
  future.
"""
struct RepeatedGame{N, T<:Real, TD<:Real}
    sg::NormalFormGame{N, T}
    delta::TD
end

# Type alias for 2 player game
"""
    RepGame2

Type representing a 2-player repeated game; alias for `RepeatedGame{2}`.
"""
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
- `delta::TD` : The common discount rate at which all players discount the
  future.

# Returns

- `::RepeatedGame` : The repeated game.
"""
RepeatedGame(p1::Player, p2::Player, delta::TD) where TD =
    RepeatedGame(NormalFormGame((p1, p2)), delta)

"""
    unpack(rpd)

Helper function that unpacks the elements of a repeated game.

# Arguments

- `rpd::RepeatedGame` : The repeated game.

# Returns

- `::Tuple{NormalFormGame, TD}` : A tuple containing the stage game
  and the delta.
"""
unpack(rpd::RepeatedGame) = (rpd.sg, rpd.delta)

"""
    _coefficient_type(rpd)

Helper function that determines the coefficient type for computations based on 
the repeated game's payoff and delta types.

# Arguments

- `rpd::RepeatedGame` : The repeated game.

# Returns

- `::Type` : Rational{BigInt} if payoffs are rational/integer and delta is rational, 
  Float64 otherwise.
"""
_coefficient_type(::RepeatedGame{N, <:Union{Rational,Integer}, <:Rational}) where {N} = Rational{BigInt}
_coefficient_type(::RepeatedGame) = Float64

# Flow utility in terms of the players actions
flow_u_1(rpd::RepGame2, a1::Int, a2::Int) =
    rpd.sg.players[1].payoff_array[a1, a2]
flow_u_2(rpd::RepGame2, a1::Int, a2::Int) =
    rpd.sg.players[2].payoff_array[a2, a1]
flow_u(rpd::RepGame2, a1::Int, a2::Int) =
    [flow_u_1(rpd, a1, a2), flow_u_2(rpd, a1, a2)]

# Computes each players best deviation given an opponent's action
best_dev_i(rpd::RepGame2, i::Int, aj::Int) =
    argmax(rpd.sg.players[i].payoff_array[:, aj])
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

Places `npts` equally spaced points along the 2 dimensional unit circle and
returns the points with x coordinates in first column and y coordinates in
second column.

# Arguments

- `npts::Int` : Number of points to be placed.

# Returns

- `pts::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the coordinates
  of the points.
"""
function unitcircle(npts::Int)
    # Want our points placed on [0, 2π]
    incr = 2π / npts
    degrees = 0.0:incr:2π

    # Points on circle
    pts = Array{Float64}(undef, npts, 2)
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

- `C::Vector{Float64}` : Vector of length `nH` containing the hyperplane
  levels.
- `H::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the subgradients.
- `Z::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the extreme
  points of the value set.
"""
function initialize_sg_hpl(nH::Int, o::Vector{Float64}, r::Float64)
    # First create unit circle
    H = unitcircle(nH)
    HT = H'

    # Choose origin and radius for big approximation
    Z = Array{Float64}(undef, 2, nH)
    for i=1:nH
        # We know that players can ever get worse than their
        # lowest punishment, so ignore anything below that
        Z[1, i] = o[1] + r*HT[1, i]
        Z[2, i] = o[2] + r*HT[2, i]
    end

    # Corresponding hyperplane levels
    C = dropdims(sum(HT .* Z, dims=1), dims=1)

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

- `C::Vector{Float64}` : Vector of length `nH` containing the hyperplane
  levels.
- `H::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the subgradients.
- `Z::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the extreme
  points of the value set.
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
- `H::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the subgradients
  used to approximate the value set, where `nH` is the number of subgradients.

# Returns

- `c::Vector{Float64}` : Vector of length `nH` used to determine which
  subgradient should be used, where `nH` is the number of subgradients.
- `A::Matrix{Float64}` : Matrix of shape `(nH+2, 2)` with nH set constraints
  and to be filled with 2 additional incentive compatibility constraints.
- `b::Vector{Float64}` : Vector of length `nH+2` to be filled with the values
  for the constraints.
"""
function initialize_LP_matrices(rpd::RepGame2, H::Matrix{Float64})
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
    b = Array{Float64}(undef, nH + 2)

    return c, A, b
end

"""
    worst_value_i(rpd, H, C, i, lp_solver=highs_optimizer_silent)

Given a constraint w ∈ W, this finds the worst possible payoff for agent i.

# Arguments

- `rpd::RepGame2` : Two player repeated game.
- `H::Matrix{Float64}` : Matrix of shape `(nH, 2)` containing the subgradients
  here `nH` is the number of subgradients.
- `C::Vector{Float64}` : The array containing the hyperplane levels.
- `i::Int` : The player of interest.
- `lp_solver` : Linear programming solver to be used internally. Pass a
  `MathOptInterface.AbstractOptimizer` type (such as `HiGHS.Optimizer`) if no
  option is needed, or a function (such as `GameTheory.highs_optimizer_silent`)
  to supply options.


# Returns

- `out::Float64` : Worst possible payoff for player i.
"""
function worst_value_i(
    rpd::RepGame2, H::Matrix{Float64},
    C::Vector{Float64}, i::Int,
    lp_solver=highs_optimizer_silent
)
    # Objective depends on which player we are minimizing
    c = zeros(2)
    c[i] = 1.0

    CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
    optimizer = MOIU.CachingOptimizer(CACHE, lp_solver())

    # Add variables
    x = MOI.add_variables(optimizer, 2)

    # Define objective function
    MOI.set(optimizer,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # Add constraints
    for i in 1:size(H,1)
        MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(H[i, :],x), 0.0),
            MOI.LessThan(C[i])
        )
    end

    # Optimize
    MOI.optimize!(optimizer)

    status = MOI.get(optimizer, MOI.TerminationStatus())

    if status == MOI.OPTIMAL
        variable_result = MOI.get(optimizer, MOI.VariablePrimal(), x)
        out = variable_result[i]
    else
        out = minimum(rpd.sg.players[i].payoff_array)
    end

    return out
end

"See `worst_value_i` for documentation"
worst_value_1(
    rpd::RepGame2,
    H::Matrix{Float64},
    C::Vector{Float64},
    lp_solver=highs_optimizer_silent
) = worst_value_i(rpd, H, C, 1, lp_solver)

"See `worst_value_i` for documentation"
worst_value_2(
    rpd::RepGame2,
    H::Matrix{Float64},
    C::Vector{Float64},
    lp_solver=highs_optimizer_silent
) = worst_value_i(rpd, H, C, 2, lp_solver)

#
# Outer Hyper Plane Approximation
#
"""
    outerapproximation(rpd; nH=32, tol=1e-8, maxiter=500, check_pure_nash=true,
                       verbose=false, nskipprint=50,
                       plib=CDDLib.Library(),
                       lp_solver=GameTheory.highs_optimizer_silent)

Approximates the set of equilibrium values for a repeated game with the outer
hyperplane approximation described by Judd, Yeltekin, Conklin (2002).

# Arguments

- `rpd::RepGame2` : Two player repeated game.
- `nH::Int` : Number of subgradients used for the approximation.
- `tol::Float64` : Tolerance in differences of set.
- `maxiter::Int` : Maximum number of iterations.
- `check_pure_nash`: Whether to perform a check about whether a pure Nash
  equilibrium exists.
- `verbose::Bool` : Whether to display updates about iterations and distance.
- `nskipprint::Int` : Number of iterations between printing information
  (assuming verbose=true).
- `plib::Polyhedra.Library`: Allows users to choose a particular package for
  the geometry computations.
  (See [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)
  docs for more info). By default, it chooses to use `CDDLib.Library()`.
- `lp_solver` : Linear programming solver to be used internally. Pass a
  `MathOptInterface.AbstractOptimizer` type (such as `HiGHS.Optimizer`) if no
  option is needed, or a function (such as `GameTheory.highs_optimizer_silent`)
  to supply options.

# Returns

- `vertices::Matrix{Float64}` : Vertices of the outer approximation of the
  value set.
"""
function outerapproximation(
        rpd::RepGame2; nH::Int=32, tol::Float64=1e-8, maxiter::Int=500,
        check_pure_nash::Bool=true, verbose::Bool=false, nskipprint::Int=50,
        plib::Polyhedra.Library=CDDLib.Library(),
        lp_solver=highs_optimizer_silent
    )

    # set up optimizer
    CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
    optimizer = MOIU.CachingOptimizer(CACHE, lp_solver())

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

    # Set iterative parameters and iterate until converged
    iter, dist = 0, 10.0
    while (iter < maxiter) & (dist > tol)
        # Compute the current worst values for each agent
        _w1 = worst_value_1(rpd, H, C, lp_solver)
        _w2 = worst_value_2(rpd, H, C, lp_solver)

        # Update all set constraints -- Copies elements 1:nH of C into b
        copyto!(b, 1, C, 1, nH)

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
            Cia = Array{Float64}(undef, nAS)
            Wia = Array{Float64}(undef, 2, nAS)
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

                MOI.empty!(optimizer)

                # Add variables
                x = MOI.add_variables(optimizer, 2)

                # Define objective function
                MOI.set(
                    optimizer,
                    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
                    )
                MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

                # Add constraints
                for i in 1:size(A,1)
                    MOI.add_constraint(optimizer,
                                       MOI.ScalarAffineFunction(
                                       MOI.ScalarAffineTerm.(A[i, :],x), 0.0),
                                       MOI.LessThan(b[i]))
                end

                # Solve corresponding linear program
                MOI.optimize!(optimizer)

                status = MOI.get(optimizer, MOI.TerminationStatus())
                if status == MOI.OPTIMAL
                    # Pull out optimal value and compute
                    w_sol = MOI.get(optimizer, MOI.VariablePrimal(), x)
                    value = (1-delta)*flow_u(rpd, a1, a2) + delta*w_sol

                    # Save hyperplane level and continuation promises
                    Cia[ia] = h1*value[1] + h2*value[2]
                    Wia[:, ia] = value
                else
                    Cia[ia] = -Inf
                end
            end

            # Action which pushes furthest in direction h_i
            astar = argmax(Cia)
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
            @warn "Maximum Iteration Reached"
        end

        # Update hyperplane levels
        copyto!(C, Cnew)
    end


    # Given the H-representation `(H, C)` of the computed polytope of
    # equilibrium payoff profiles, we obtain its V-representation `vertices`
    # using Polyhedra.jl (it uses `plib` which was chosen for computations)
    p = polyhedron(hrep(H, C), plib)
    vr = vrep(p)
    pts = points(vr)  # Vector of Vectors

    # Reduce the number of vertices by rounding points to the tolerance
    tol_int = round(Int, abs(log10(tol))) - 1

    # Find vertices that are unique within tolerance level
    vertices = Matrix{Float64}(undef, (length(pts), 2))
    for (i, pt) in enumerate(pts)
        vertices[i, :] = round.(pt, digits=tol_int)
    end
    vertices = unique(vertices, dims=1)

    return vertices
end

"""
    _find_ic_boundary_points_optimization(p, IC1, IC2, S)

Find intersection points on IC boundaries using optimization instead of 
polyhedron intersection. This implements optimization 2 from PR #65 comment 5.

For a given polyhedron p and IC constraints IC1, IC2, find points where:
- v[1] = IC1 (vertical IC boundary)  
- v[2] = IC2 (horizontal IC boundary)

Uses linear programming to find extreme points on these boundaries.
"""
function _find_ic_boundary_points_optimization(p::Polyhedron, IC1::T, IC2::T, S::Type) where T<:Real
    boundary_points = Vector{Vector{S}}()
    
    # Get the H-representation of the polyhedron
    h = hrep(p)
    if isempty(h)
        return boundary_points
    end
    
    try
        # Extract constraint matrix A and vector b from Ax <= b
        inequalities = collect(halfspaces(h))
        if isempty(inequalities)
            return boundary_points
        end
        
        A = zeros(S, length(inequalities), 2)
        b = zeros(S, length(inequalities))
        
        for (i, ineq) in enumerate(inequalities)
            # inequality is of form a'x <= β
            A[i, :] = ineq.a
            b[i] = ineq.β
        end
        
        # Setup optimizer
        optimizer = HiGHS.Optimizer()
        MOI.set(optimizer, MOI.Silent(), true)
        
        # Find points on v[1] = IC1 boundary (maximize and minimize v[2])
        for objective_sense in [MOI.MAX_SENSE, MOI.MIN_SENSE]
            model = MOI.instantiate(optimizer)
            
            # Variables: v = [v1, v2]
            v = MOI.add_variables(model, 2)
            
            # Constraint: v[1] = IC1
            MOI.add_constraint(model, v[1], MOI.EqualTo(S(IC1)))
            
            # Constraints from polyhedron: A*v <= b
            for i in 1:size(A, 1)
                MOI.add_constraint(model, 
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[i, :], v), zero(S)),
                    MOI.LessThan(b[i]))
            end
            
            # Objective: maximize/minimize v[2] 
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{S}}(),
                    MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(one(S), v[2])], zero(S)))
            MOI.set(model, MOI.ObjectiveSense(), objective_sense)
            
            MOI.optimize!(model)
            
            if MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
                solution = MOI.get(model, MOI.VariablePrimal(), v)
                # Check if this point is approximately on the IC1 boundary
                if abs(solution[1] - IC1) < 1e-10
                    push!(boundary_points, S[solution[1], solution[2]])
                end
            end
        end
        
        # Find points on v[2] = IC2 boundary (maximize and minimize v[1])
        for objective_sense in [MOI.MAX_SENSE, MOI.MIN_SENSE]
            model = MOI.instantiate(optimizer)
            
            # Variables: v = [v1, v2]
            v = MOI.add_variables(model, 2)
            
            # Constraint: v[2] = IC2
            MOI.add_constraint(model, v[2], MOI.EqualTo(S(IC2)))
            
            # Constraints from polyhedron: A*v <= b
            for i in 1:size(A, 1)
                MOI.add_constraint(model, 
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[i, :], v), zero(S)),
                    MOI.LessThan(b[i]))
            end
            
            # Objective: maximize/minimize v[1]
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{S}}(),
                    MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(one(S), v[1])], zero(S)))
            MOI.set(model, MOI.ObjectiveSense(), objective_sense)
            
            MOI.optimize!(model)
            
            if MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
                solution = MOI.get(model, MOI.VariablePrimal(), v)
                # Check if this point is approximately on the IC2 boundary
                if abs(solution[2] - IC2) < 1e-10
                    push!(boundary_points, S[solution[1], solution[2]])
                end
            end
        end
    catch e
        # Fall back to empty result if optimization fails
        return boundary_points
    end
    
    # Remove duplicates
    unique_points = Vector{Vector{S}}()
    for pt in boundary_points
        is_duplicate = false
        for existing in unique_points
            if norm(pt - existing) < 1e-10
                is_duplicate = true
                break
            end
        end
        if !is_duplicate
            push!(unique_points, pt)
        end
    end
    
    return unique_points
end

"""
    _add_point_incremental!(p_current, new_point, lib, S)

Add a point to the current polyhedron incrementally. This implements 
optimization 1 from PR #65 comment 5.

Returns updated polyhedron if point was added, or original polyhedron if point was redundant.
"""
function _add_point_incremental!(p_current::Polyhedron, new_point::AbstractVector, lib, S::Type)
    # Check if point is already in the polyhedron
    h = hrep(p_current)
    if new_point in h
        return p_current  # Point is redundant, don't add
    end
    
    # Add the point to V-representation
    v_current = vrep(p_current)
    pts_current = points(v_current)
    
    # Create new points array with the additional point
    new_points = [collect(pt) for pt in pts_current]
    push!(new_points, collect(new_point))
    
    # Create new polyhedron
    p_new = polyhedron(vrep(new_points), lib)
    removevredundancy!(p_new)
    
    return p_new
end

"""
    AS(rpd; maxiter=1000, plib=default_library(2, Float64), tol=1e-5, u=nothing, verbose=false, use_optimization=false, incremental_redundancy=false)

Using AS algorithm to compute the set of payoff pairs of all pure-strategy
subgame-perfect equilibria with public randomization for any repeated
two-player games with perfect monitoring and discounting, following
Abreu and Sannikov (2014).

# Arguments

- `rpd::RepeatedGame{2, T, TD}` : Two player repeated game with T<:Real, TD<:Real.
- `maxiter::Integer` : Maximum number of iterations.
- `plib`: Allows users to choose a particular package for the geometry
  computations.
  (See [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)
  docs for more info). By default, it chooses to use SimplePolyhedraLibrary.
- `tol::Float64` : Tolerance in differences of set.
- `u` : The punishment payoff pair if any player deviates. In default,
  we use minimax payoff pair. If there is better guess, you can specify it
  by passing a `Vector` with length 2.
- `verbose::Bool` : If true, print convergence information. Defaults to false.
- `use_optimization::Bool` : If true, use optimization-based IC boundary intersection finding (PR #65 suggestion 2). Defaults to false.
- `incremental_redundancy::Bool` : If true, use incremental redundancy removal (PR #65 suggestion 1). Defaults to false.

# Returns

- `::Matrix{S}` : Vertices of the set of payoff pairs, where S is determined by _coefficient_type(rpd).
"""
function AS(rpd::RepeatedGame{2,T,TD}; maxiter::Integer=1000,
            plib=default_library(2, Float64), tol::Float64=1e-5,
            u::Union{AbstractVector, Nothing}=nothing, verbose::Bool=false,
            use_optimization::Bool=false, incremental_redundancy::Bool=false) where {T,TD}

    S = _coefficient_type(rpd)
    lib = similar_library(plib, 2, S)

    # Initialize W0 with each entries of payoff bimatrix
    v_old = _payoff_points(S, rpd.sg)

    if isnothing(u)
        u = S[minimum(rpd.sg.players[1].payoff_array),
                    minimum(rpd.sg.players[2].payoff_array)]
    else
        u = convert(Vector{S}, u)
    end

    # create VRepresentation and Polyhedron and get rid of redundant vertices
    p = polyhedron(vrep(v_old), lib)
    removevredundancy!(p)
    H = hrep(p)

    # calculate the best deviation gains
    # normalize with (1-delta)/delta
    best_dev_gains1, best_dev_gains2 = (one(S)-rpd.delta)/rpd.delta .* _best_dev_gains(rpd.sg)

    for iter = 1:maxiter

        if incremental_redundancy
            # Optimization 1: Use incremental polyhedron construction
            p_current = deepcopy(p)
        else
            v_new = S[] # to store new vertices
            # Use sizehint! for better performance as suggested in PR #65
            sizehint!(v_new, 8 * prod(rpd.sg.nums_actions))
        end
        
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
                        if incremental_redundancy
                            p_current = _add_point_incremental!(p_current, S[payoff1, payoff2], lib, S)
                        else
                            push!(v_new, payoff1, payoff2)
                        end
                    end
                end

                # find out the intersections of polyhedron and IC boundaries
                if use_optimization
                    # Optimization 2: Use optimization-based boundary finding
                    try
                        boundary_points = _find_ic_boundary_points_optimization(p, IC1, IC2, S)
                        if !isempty(boundary_points)
                            for pt in boundary_points
                                new_point = rpd.delta * pt + (one(S) - rpd.delta) * S[payoff1, payoff2]
                                if incremental_redundancy
                                    p_current = _add_point_incremental!(p_current, new_point, lib, S)
                                else
                                    push!(v_new, new_point...)
                                end
                            end
                        else
                            # Fallback to original method if optimization fails to find points
                            p_IC = polyhedron(hrep(-Matrix{S}(I, 2, 2), -S[IC1, IC2]), lib)
                            p_inter = intersect(p_IC, p)
                            Vmat = MixedMatVRep(vrep(p_inter)).V
                            for i in 1:size(Vmat, 1)
                                if Vmat[i, 1] ≈ IC1 || Vmat[i, 2] ≈ IC2
                                    new_point = (rpd.delta * Vmat[i, :] +
                                                  (one(S) - rpd.delta) * S[payoff1, payoff2])
                                    if incremental_redundancy
                                        p_current = _add_point_incremental!(p_current, new_point, lib, S)
                                    else
                                        push!(v_new, new_point...)
                                    end
                                end
                            end
                        end
                    catch
                        # Fallback to original method if optimization completely fails
                        p_IC = polyhedron(hrep(-Matrix{S}(I, 2, 2), -S[IC1, IC2]), lib)
                        p_inter = intersect(p_IC, p)
                        Vmat = MixedMatVRep(vrep(p_inter)).V
                        for i in 1:size(Vmat, 1)
                            if Vmat[i, 1] ≈ IC1 || Vmat[i, 2] ≈ IC2
                                new_point = (rpd.delta * Vmat[i, :] +
                                              (one(S) - rpd.delta) * S[payoff1, payoff2])
                                if incremental_redundancy
                                    p_current = _add_point_incremental!(p_current, new_point, lib, S)
                                else
                                    push!(v_new, new_point...)
                                end
                            end
                        end
                    end
                else
                    # Original method: construct p_IC and intersect
                    p_IC = polyhedron(hrep(-Matrix{S}(I, 2, 2), -S[IC1, IC2]), lib)
                    p_inter = intersect(p_IC, p)
                    Vmat = MixedMatVRep(vrep(p_inter)).V
                    for i in 1:size(Vmat, 1)
                        if Vmat[i, 1] ≈ IC1 || Vmat[i, 2] ≈ IC2
                            new_point = (rpd.delta * Vmat[i, :] +
                                          (one(S) - rpd.delta) * S[payoff1, payoff2])
                            if incremental_redundancy
                                p_current = _add_point_incremental!(p_current, new_point, lib, S)
                            else
                                push!(v_new, new_point...)
                            end
                        end
                    end
                end
            end
        end

        if incremental_redundancy
            # Use the incrementally built polyhedron
            p = p_current
            v_dedup = MixedMatVRep(vrep(p)).V
        else
            v_new = reshape(v_new, 2, :)'

            # get rid of redundant points
            p = polyhedron(vrep(v_new), lib)
            removevredundancy!(p)

            # check if it's converged
            # Use deduplicated vertices for convergence check
            v_dedup = MixedMatVRep(vrep(p)).V
        end

        # check if it's converged
        # first check if the numbers of vertices are the same
        if size(v_dedup) == size(v_old)
            # then check the euclidean distance
            if norm(v_dedup-v_old) < tol
                verbose && println("converged in $(iter) iterations")
                break
            end
        end

        # check if maxiter is reached
        if iter == maxiter
            @warn "Maximum Iteration Reached"
        end

        v_old = v_dedup
        H = hrep(p)

        # step 2
        # update u
        if incremental_redundancy
            # Extract vertices from the polyhedron for minimum calculation
            v_matrix = MixedMatVRep(vrep(p)).V
            u_ = [minimum(v_matrix[:, 1]),
                  minimum(v_matrix[:, 2])]
        else
            u_ = [minimum(v_new[:, 1]),
                  minimum(v_new[:, 2])]
        end
        if any(u_ .> u)
            u = u_
        end
    end

    # Return matrix with coefficient type S
    vr = vrep(p)
    pts = points(vr)
    vertices = Matrix{S}(undef, (length(pts), 2))
    for (i, pt) in enumerate(pts)
        vertices[i, :] = S.(pt)
    end

    return vertices
end

"""
    uniquetolrows(V, tol)

Remove near-duplicate rows from matrix V using tolerance-based deduplication.

# Arguments

- `V::AbstractMatrix{T}` : Input matrix where T<:Real.
- `tol::Real` : Tolerance for considering points as duplicates.

# Returns

- `::Matrix{T}` : Matrix with duplicate rows removed within tolerance.
"""
function uniquetolrows(V::AbstractMatrix{T}, tol::Real) where T
    digits = max(0, floor(Int, -log10(tol)))
    Vr = round.(V; digits=digits)
    return unique(Vr; dims=1)
end
"""
    _payoff_points(::Type{T}, g)

Return a matrix with each row being a payoff pair point in the two dimensional
space.

# Arguments

- `g::NormalFormGame{2}` : Two-player NormalFormGame.

# Returns

- `v::Matrix{T}` : Matrix with size n by 2, where n is the number of
  action profiles. Each row corresponds to one payoff pair.
"""
function _payoff_points(::Type{T}, g::NormalFormGame{2}) where T

    nums_action_profiles = prod(g.nums_actions)
    v = Matrix{T}(undef, nums_action_profiles, 2)
    v[:, 1] = reshape(g.players[1].payoff_array, nums_action_profiles)
    v[:, 2] = reshape(g.players[2].payoff_array', nums_action_profiles)

    return v
end

"""
    _best_dev_gains(g)

Calculate the payoff gains from deviating from the current action to
the best response for each player.

# Arguments

- `g::NormalFormGame{2, T}` : Two-player NormalFormGame.

# Returns

- `::Tuple{Matrix{T}, Matrix{T}}` : Tuple of best deviating gain matrices
  for two players. For example, for the first matrix `best_dev_gains1`,
  `best_dev_gains1[i, j]` is the payoff gain for player 1 for deviating
  to the best response from ith action given player 2 choosing jth action.
"""
function _best_dev_gains(g::NormalFormGame{2, T}) where T

    best_dev_gains1 = (maximum(g.players[1].payoff_array; dims=1)
                       .- g.players[1].payoff_array)
    best_dev_gains2 = (maximum(g.players[2].payoff_array; dims=1)
                       .- g.players[2].payoff_array)

    return best_dev_gains1, best_dev_gains2
end
