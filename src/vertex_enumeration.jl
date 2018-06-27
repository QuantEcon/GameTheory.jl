import Polyhedra: hyperplane

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

function _vertex_enumeration_producer{T}(c::Channel,
                                         g::NormalFormGame{2, T},
                                         plib)

    n, m = size(g.players[1].payoff_array)

    # create Representation for player 1
    H1, V1, H2, V2 = construction_BRP(g, plib)

    ZERO_LABELING_BITS = (1 << (n+m)) - (1 << m)
    COMPLETE_LABELING_BITS = 1 << (n+m) - 1

    for v1 in points(V1)
        labelings_bits1 = labelings_bits(v1, H1)
        if labelings_bits1 == ZERO_LABELING_BITS
            continue
        end
        for v2 in points(V2)
            labelings_bits2 = labelings_bits(v2, H2)
            if xor(labelings_bits1, labelings_bits2) == COMPLETE_LABELING_BITS
                put!(c, (_get_mixed_action(v1),
                         _get_mixed_action(v2)))
            end
        end
    end

end

function construction_BRP{T}(g::NormalFormGame{2, T}, plib)

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
    V1 = vrep(p1)

    # create Representation for player 2
    D = Matrix{T}(n+m, m)
    D[1:m, :] = -eye(T, m)
    D[m+1:end, :] = g.players[1].payoff_array
    b2 = Vector{T}(n+m)
    b2[1:m] = zero(T)
    b2[m+1:end] = one(T)
    H2 = hrep(D, b2)
    p2 = polyhedron(H2, plib)
    V2 = vrep(p2)

    return H1, V1, H2, V2
end

function labelings_bits(v::VRepElement, p::HRep)
    b = 0
    for (i, h) in enumerate(halfspaces(p))
        if v in hyperplane(h)
            b += 1 << (i-1)
        end
    end
    return b
end

function _get_mixed_action(a::Vector)
    return a ./ sum(a)
end

