#=
Utility functions used in Games.jl
=#

"""
    clp_optimizer_silent()

Function that returns a `Clp.Optimizer` instance with `LogLevel` set to `0`.
"""
function clp_optimizer_silent()
    optimizer = Clp.Optimizer()
    MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end
