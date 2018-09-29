using Games
using Compat
using Compat.Test

include("test_pure_nash.jl")
include("test_repeated_game.jl")
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")

include("generators/runtests.jl")
