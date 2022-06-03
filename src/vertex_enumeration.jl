function nonnegativeorthant_hrep(dim::Int)
    h = Vector{HalfSpace{Float64, Vector{Float64}}}()
    for i in 1:dim
        e_i = zeros(dim)
        e_i[i] = -1.
        push!(h, HalfSpace(e_i, 0.))
    end
    H = h[1]
    for i in 2:lastindex(h)
        H = H ∩ h[i]
    end
    H
end
nonnegativeorthant(dim::Int) = polyhedron(nonnegative_hrep(dim))


function br_envelope_hrep(player::Player)
    A = player.payoff_array
    h = Vector{HalfSpace{Float64, Vector{Float64}}}()
    for i in 1:num_actions(player)
        A_i = A[i,:]
        push!(h, HalfSpace(A_i, 1.))
    end
    H = h[1]
    for i in 2:lastindex(h)
        H = H ∩ h[i]
    end
    H
end


function bestresponsepolyhedra(G::NormalFormGame)
    nnorthantP = nonnegativeorthant_hrep(num_actions(G.players[1]))  # x_i ≥ 0 for all i=1,…,m
    nnorthantQ = nonnegativeorthant_hrep(num_actions(G.players[2]))  # x_i ≥ 0 for all i=m+1,…,m+n
    brenvelopeP = br_envelope_hrep(G.players[2])  # B'x ≤ 1
    brenvelopeQ = br_envelope_hrep(G.players[1])  # Ax ≤ 1
    P = brenvelopeP ∩ nnorthantP  # best response polyhedron for Player 1
    Q = nnorthantQ ∩ brenvelopeQ  # best response polyhedron for Player 2
    polyhedron(P), polyhedron(Q)
end


function hlabels(P::HRepresentation)
    Phindices = Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}[]
    for pi in eachindex(halfspaces(P))
        push!(Phindices, pi)
    end
    Phindices
end
hlabels(P::DefaultPolyhedron) = hlabels(hrep(P))
hlabels(P::VRepresentation) = hlabels(doubledescription(P))


function vlabels(P::VRepresentation)
    Pvindices = Polyhedra.Index{Float64, Vector{Float64}}[]
    for pi in eachindex(points(P))
        push!(Pvindices, pi)
    end
    Pvindices
end
vlabels(P::DefaultPolyhedron) = vlabels(vrep(P))
vlabels(P::HRepresentation) = vlabels(doubledescription(P))


label_to_integer(idx::Polyhedra.Index{T, S}) where {T, S} = idx.value
# integer_to_vlabel(n::Int) = Polyhedra.Index{Float64, Vector{Float64}}(n)
# integer_to_hlabel(n::Int) = Polyhedra.Index{Float64, HalfSpace{Float64, Vector{Float64}}}(n)


function labelmap(P::DefaultPolyhedron)
    labelmaps = []
    vpoints = [x for x in points(P)]
    for pi in eachindex(points(P))
        push!(labelmaps, (vpoints[label_to_integer(pi)], incidenthalfspaceindices(P, pi)))
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


function LabeledBimatrixGame(G::NormalFormGame)
    @assert num_players(G) == 2
    P, Q = bestresponsepolyhedra(G)
    LabeledBimatrixGame(G, LabeledPolyhedron(P), LabeledPolyhedron(Q))
end


# unlabel(LP::LabeledPolyhedron) = LP.polyhedron
# unlabel(P::T) where T <: Union{Polyhedron, HRepresentation, VRepresentation} = P


get_num_hlabels(LP::LabeledPolyhedron, point) = length(LP.labelmap[point])


function is_nondegenerate(B::LabeledBimatrixGame)
    m = num_actions(B.game.players[1])
    n = num_actions(B.game.players[2])
    all([(get_num_hlabels(B.P, point) ≤ m) for point in B.P.points]) && all([(get_num_hlabels(B.Q, idx) ≤ n) for idx in B.Q.points])
end


function is_nondegenerate(G::NormalFormGame)
    if num_players(G) == 2
        return is_nondegenerate(LabeledBimatrixGame(G))
    else
        error("Not a bimatrix game")
    end
end


#function droplabel!(LP::LabeledPolyhedron, point, label)
#    LP.labelmap[point]
#    if !(label in LP.labelmap[point])
#        @warn("Label not found.")
#        return nothing
#    end
#    deleteat!(LP.labelmap[point], findall(x->x==label, LP.labelmap[point])[1])  # note label is unique so findall() finds precisely 1 index
#end


#function droplabel(LP::LabeledPolyhedron, point, label)
#    LP.labelmap[point]
#    if !(label in LP.labelmap[point])
#        @warn("Label not found.")
#        return LP.labelmap[point]
#    end
#    deleteat(LP.labelmap[point], findall(x->x==label, LP.labelmap[point])[1])  # note label is unique so findall() finds precisely 1 index
#end


function dropvertex_pure!(LP::LabeledPolyhedron, point)
    delete!(LP.labelmap, point)
    deleteat!(LP.points, findall(x->x==point, LP.points)[1])
end

function dropvertex!(LP::LabeledPolyhedron, point)
    if length(findall(x -> x==point, LP.points)) == 0
        dropvertex!(LP::LabeledPolyhedron, point, 1e-10)
    else
        dropvertex_pure!(LP, point)
    end
end


function dropvertex!(LP::LabeledPolyhedron, point, tol)
    local_point = LP.points[findall(x->norm(x-point)<tol, LP.points)[1]]
    local_point
    dropvertex_pure!(LP, local_point)
end


"""
    vertex_enumeration(G::NormalFormGame)

Finds all Nash equilibria via the equilibria by vertex enumeration algorithm (Algorithm 3.5 in von Stengel (2007).)
"""
function vertex_enumeration(g::NormalFormGame)
    B = LabeledBimatrixGame(g)
    if !is_nondegenerate(B)
        @error("The vertex enumeration algorithm will not yield a solution for degenerate games.")
        return nothing
    end
    m, n = num_actions.(g.players)
    nash = []
    dropvertex!(B.P, zeros(m))
    dropvertex!(B.Q, zeros(n))
    if length(B.P.points) ≤ 1 || length(B.Q.points) ≤ 1
        @warn("At least one best response polyhedron is collapsed to the origin. Applying positive affine transformations.")
        player1_transform = Player(G.players[1].payoff_array .- minimum(G.players[1].payoff_array))
        player2_transform = Player(G.players[2].payoff_array .- minimum(G.players[2].payoff_array))
        B = LabeledBimatrixGame(NormalFormGame(player1_transform, player2_transform))
        dropvertex!(B.P, zeros(m))
        dropvertex!(B.Q, zeros(n))
    end
    for x in B.P.points
        for y in B.Q.points
            if sort(label_to_integer.(vcat(B.P.labelmap[x], B.Q.labelmap[y])))== Vector(1:m+n) 
                # i.e. (x, y) completely labeled, x ∈ P - {0}, y ∈ Q - {0} 
                push!(nash, (x./(ones(1,m)*x),y./(ones(1,n)*y)))
            end
        end
    end
    nash
end