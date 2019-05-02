using Games
using Test

include("test_pure_nash.jl")
include("test_repeated_game.jl")
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")
include("test_fictplay.jl")
include("test_localint.jl")
include("test_brd.jl")
include("test_logitdyn.jl")

include("generators/runtests.jl")
