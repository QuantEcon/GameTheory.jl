using Games
using Compat
using Base.Test

include("test_pure_nash.jl")
if VERSION >= v"0.7-"
	@test_skip include("test_repeated_game.jl")
else
	include("test_repeated_game.jl")
end
include("test_normal_form_game.jl")
include("test_random.jl")
include("test_support_enumeration.jl")

include("generators/runtests.jl")
