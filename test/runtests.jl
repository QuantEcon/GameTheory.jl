using Games

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("test_normal_form_game.jl")
include("test_pure_nash.jl")

