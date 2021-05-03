function vertex_enumeration(g::NormalFormGame{2};
                            plib=SimplePolyhedraLibrary{Float64}())

    c = Channel(0)
    task = vertex_enumeration_task(c, g, plib)
    bind(c, task)
    schedule(task)
    NEs = Tuple{Vector{Real}, Vector{Real}}[NE for NE in c]

    return NEs

end

function vertex_enumeration_task(c::Channel,
                                 g::NormalFormGame{2},
                                 plib)

    task = Task(
        () -> _vertex_enumeration_producer(c, g, plib)
    )

    return task

end

function _vertex_enumeration_producer(c::Channel,
                                      g::NormalFormGame{2, T},
                                      plib)

    n, m = size(g.players[1].payoff_array)

    # create Representation for player 1
    p1, p2 = construction_BRP(g, plib)
    V1 = points(p1)
    simplex1 = []
    for pidx in eachindex(points(p1))
        push!(simplex1, [idx.value for idx in incidenthalfspaceindices(p1, pidx)])
    end
    
    V2 = points(p2)
    simplex2 = []
    for pidx in eachindex(points(p2))
        push!(simplex2, [idx.value for idx in incidenthalfspaceindices(p2, pidx)])
    end

    ZERO_LABELING_BITS = (1 << (n+m)) - (1 << m)
    COMPLETE_LABELING_BITS = 1 << (n+m) - 1
    
    for (i, v1) in enumerate(V1)
        labelings_bits1 = labelings_bits(simplex1[i])
        if labelings_bits1 == ZERO_LABELING_BITS
            continue
        end
        for (j, v2) in enumerate(V2)
            labelings_bits2 = labelings_bits(simplex2[j])
            if xor(labelings_bits1, labelings_bits2) == COMPLETE_LABELING_BITS
                put!(c, (_get_mixed_action(v1),
                         _get_mixed_action(v2)))
            end
        end
    end

end

function construction_BRP(g::NormalFormGame{2, T}, plib)

    n, m = size(g.players[1].payoff_array)

    # create Representation for player 1
    C = Matrix{T}(n+m, n)
    C[1:m, :] = g.players[2].payoff_array
    C[m+1:end, :] = -eye(T, n)
    b1 = Vector{T}(n+m)
    b1[1:m] = one(T)
    b1[m+1:end] = zero(T)
    H1 = hrep(C, b1)
    p1 = polyhedron(H1, plib)

    # create Representation for player 2
    D = Matrix{T}(n+m, m)
    D[1:m, :] = -eye(T, m)
    D[m+1:end, :] = g.players[1].payoff_array
    b2 = Vector{T}(n+m)
    b2[1:m] = zero(T)
    b2[m+1:end] = one(T)
    H2 = hrep(D, b2)
    p2 = polyhedron(H2, plib)

    return p1, p2
end

function labelings_bits(inchalfindices::Vector{T}) where T <: Integer
    b = 0
    for i in inchalfindices
        b += 1 << (i-1)
    end
    return b
end

function _get_mixed_action(a::Vector{T})
    return a ./ sum(a)
end
