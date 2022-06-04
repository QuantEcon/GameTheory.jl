using GameTheory
using Test

include("test_pure_nash.jl")
include("test_repeated_game.jl")
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")
include("test_lrsnash.jl")
include("test_vertex_enumeration.jl")


include("generators/runtests.jl")
