# ------------------------------------- #
# Testing the code quality with Aqua.jl #
# ------------------------------------- #

using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        GameTheory;
        ambiguities=(recursive=true,),
        # Signatures such as `actions::NTuple{N,Vector{T}}` with `N` the number
        # of players leave `T` unbound only for `N = 0`, which cannot occur;
        # the check cannot know that and reports them
        unbound_args=false,
        # TreeViews, a dependency of HomotopyContinuation, ships without a
        # Project.toml, which this check cannot handle
        persistent_tasks=false,
    )
end
