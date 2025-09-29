using GameTheory
using Test

include("util.jl")

include("test_pure_nash.jl")
include("test_repeated_game.jl")
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")
include("test_lemke_howson.jl")
include("test_lrsnash.jl")
include("test_homotopy_continuation.jl")
include("test_fictplay.jl")
include("test_localint.jl")
include("test_brd.jl")
include("test_logitdyn.jl")

include("generators/runtests.jl")
