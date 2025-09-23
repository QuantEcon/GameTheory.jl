#=
Utility functions used in GameTheory.jl
=#

"""
    clp_optimizer_silent()

Function that returns a `Clp.Optimizer` instance in silent mode.
"""
function clp_optimizer_silent()
    optimizer = Clp.Optimizer()
    MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end
