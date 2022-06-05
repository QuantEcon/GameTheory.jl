using GameTheory
using Test

function NEs_approx_equal(NEs1::Vector{NTuple{2,Vector{T1}}},
    NEs2::Vector{NTuple{2,Vector{T2}}}) where {T1,T2}
        @test length(NEs1) == length(NEs2)
        @test T1 == T2
        for (actions1, actions2) in zip(NEs1, NEs2)
        for (action1, action2) in zip(actions1, actions2)
            @test action1 ≈ action2
        end
    end
end

include("test_pure_nash.jl")
include("test_repeated_game.jl")
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")
include("test_lrsnash.jl")
include("test_vertex_enumeration.jl")


include("generators/runtests.jl")
