#=
Entry point for the generators benchmarks

The bimatrix_generators subgroup is excluded from the default suite in
benchmarks.jl (see README.md). This script exposes it as its own `SUITE`
so that it can be run on its own, standalone:

    julia --project=benchmark benchmark/generators.jl

or through PkgBenchmark, e.g. for cross-commit comparisons:

    using PkgBenchmark
    jud = judge("GameTheory", "<target>", "<baseline>";
                script="benchmark/generators.jl")
=#
using BenchmarkTools

const SUITE = BenchmarkGroup()

SUITE["bimatrix_generators"] = include("bimatrix_generators.jl")

#= Standalone execution =#

if abspath(PROGRAM_FILE) == @__FILE__
    tune!(SUITE)
    results = run(SUITE; verbose=true)
    show(IOContext(stdout, :compact => false), MIME"text/plain"(), results)
    println()
end
