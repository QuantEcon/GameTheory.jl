#=
Utility functions used in GameTheory.jl
=#

"""
    highs_optimizer_silent()

Function that returns a `HiGHS.Optimizer` instance in silent mode.
"""
function highs_optimizer_silent()
    optimizer = HiGHS.Optimizer()
    MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end
