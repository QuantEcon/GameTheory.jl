function nonnegativeorthant_hrep(dim::Int)
    h = Vector{HalfSpace{Rational{Int64}, Vector{Rational{Int64}}}}()
    for i in 1:dim
        e_i = zeros(dim)
        e_i[i] = -1//1
        push!(h, HalfSpace(e_i, 0//1))
    end
    H = h[1]
    for i in 2:lastindex(h)
        H = H ∩ h[i]
    end
    H
end


function br_envelope_hrep(player::Player)
    A = player.payoff_array
    h = Vector{HalfSpace}()
    for i in 1:num_actions(player)
        A_i = A[i,:]
        push!(h, HalfSpace(A_i, 1))
    end
    H = h[1]
    for i in 2:lastindex(h)
        H = H ∩ h[i]
    end
    H
end


function bestresponsepolyhedra(g::NormalFormGame; plib::Polyhedra.Library =
    default_library(2, Float64))
    nnorthantP = nonnegativeorthant_hrep(num_actions(g.players[1]))  
    # x_i ≥ 0 for all i=1,…,m
    nnorthantQ = nonnegativeorthant_hrep(num_actions(g.players[2]))  
    # x_i ≥ 0 for all i=m+1,…,m+n
    brenvelopeP = br_envelope_hrep(g.players[2])  # B'x ≤ 1
    brenvelopeQ = br_envelope_hrep(g.players[1])  # Ax ≤ 1
    P = brenvelopeP ∩ nnorthantP  # best response polyhedron for Player 1
    Q = nnorthantQ ∩ brenvelopeQ  # best response polyhedron for Player 2
    polyhedron(P, plib), polyhedron(Q, plib)
end


function hlabels(P::HRepresentation)
    Phindices = Polyhedra.Index[]
    for pi in eachindex(halfspaces(P))
        push!(Phindices, pi)
    end
    Phindices
end
hlabels(P::DefaultPolyhedron) = hlabels(hrep(P))


label_to_integer(idx::Polyhedra.Index{T, S}) where {T, S} = idx.value


function labelmap(P::DefaultPolyhedron)
    labelmaps = []
    vpoints = [x for x in points(P)]
    for pi in eachindex(points(P))
        push!(labelmaps, (vpoints[label_to_integer(pi)], 
        incidenthalfspaceindices(P, pi)))
    end
    Dict(labelmaps)
end


mutable struct LabeledPolyhedron{S<:Polyhedron, T<:Vector, U<:Vector, V<:Dict}
    polyhedron::S
    points::T
    hlabels::U
    labelmap::V
end


function LabeledPolyhedron(P::DefaultPolyhedron)
    vpoints = [x for x in points(P)]
    LabeledPolyhedron(P, vpoints, hlabels(P), labelmap(P))
end


struct LabeledBimatrixGame
    game::NormalFormGame
    P::LabeledPolyhedron
    Q::LabeledPolyhedron
end


function LabeledBimatrixGame(g::NormalFormGame{2}; plib::Polyhedra.Library =
    default_library(2, Float64))
    P, Q = bestresponsepolyhedra(g; plib = plib)
    LabeledBimatrixGame(g, LabeledPolyhedron(P), LabeledPolyhedron(Q))
end


get_num_hlabels(LP::LabeledPolyhedron, point) = length(LP.labelmap[point])


function is_nondegenerate(b::LabeledBimatrixGame)
    m = num_actions(b.game.players[1])
    n = num_actions(b.game.players[2])
    all([(get_num_hlabels(b.P, point) ≤ m) for point in b.P.points]) && all([
        (get_num_hlabels(b.Q, idx) ≤ n) for idx in b.Q.points])
end


function is_nondegenerate(g::NormalFormGame)
    if num_players(g) == 2
        return is_nondegenerate(LabeledBimatrixGame(g))
    else
        error("Not a bimatrix game.")
    end
end


function dropvertex_pure!(LP::LabeledPolyhedron, point)
    delete!(LP.labelmap, point)
    deleteat!(LP.points, findall(x -> x == point, LP.points)[1])
end


function dropvertex!(LP::LabeledPolyhedron, point)
    if length(findall(x -> x == point, LP.points)) == 0
        dropvertex!(LP::LabeledPolyhedron, point, 1e-10)
    else
        dropvertex_pure!(LP, point)
    end
end


function dropvertex!(LP::LabeledPolyhedron, point, tol)  # handle precision 
    # errors up to tolerance
    local_point = LP.points[findall(x -> norm(x - point) < tol, LP.points)[1]]
    local_point
    dropvertex_pure!(LP, local_point)
end


"""
    vertex_enumeration(g::NormalFormGame{2, T})

Finds all Nash equilibria of a non-degenerate bimatrix game `g` via the vertex 
enumeration algorithm (Algorithm 3.5 in von Stengel (2007).)

# References
- B. von Stengel, "Equilibria Computation for Two-Player Games in Strategic and 
Extensive Form." 
In N. Nisan, T. Roughgarden, E. Tardos and V. V. Vazirani (eds.), Algorithmic 
Game Theory, 2007. 
"""
function vertex_enumeration(g::NormalFormGame{2, T}; plib::Polyhedra.Library = 
    default_library(2, Float64)) where T<:Real
    if T == Rational{Int64} && plib == default_library(2, Float64)
        plib = default_library(2, Rational{Int64}) # change default if payoffs 
        # are rational
    end
    if !(all(g.players[1].payoff_array .≥ 0) && all(g.players[2].payoff_array
         .≥ 0))
        player1_transform = Player(g.players[1].payoff_array 
        .- minimum(g.players[1].payoff_array))
        player2_transform = Player(g.players[2].payoff_array 
        .- minimum(g.players[2].payoff_array))
        g = NormalFormGame(player1_transform, player2_transform)  # positive 
        # affine transformation
    end
    b = LabeledBimatrixGame(g; plib = plib)
    if !is_nondegenerate(b)
        error("The vertex enumeration algorithm will not yield a solution "
        * "for degenerate games.")
    end
    m, n = num_actions.(g.players)
    NEs = NTuple{2, Vector{T}}[]
    dropvertex!(b.P, zeros(m))
    dropvertex!(b.Q, zeros(n))

    for x in b.P.points
        for y in b.Q.points
            if sort(label_to_integer.(vcat(b.P.labelmap[x], 
                b.Q.labelmap[y]))) == Vector(1:m+n) 
                # i.e. (x, y) completely labeled, x ∈ P - {0}, y ∈ Q - {0} 
                push!(NEs, (x./(ones(Rational, 1,m)*x),y./(ones(Rational, 1,n)*y)))
            end
        end
    end
    NEs
end
function vertex_enumeration(g::NormalFormGame{2, Int}; plib::Polyhedra.Library = 
    default_library(2, Float64))
    vertex_enumeration(NormalFormGame{2, Float64}(g); plib = plib)
end