#=
PkgBenchmark entry point for the generators benchmarks

The bimatrix_generators subgroup is excluded from the default suite in
benchmarks.jl (see README.md). This script exposes it as its own `SUITE`
so that PkgBenchmark can run it, e.g. for cross-commit comparisons:

    using PkgBenchmark
    jud = judge("GameTheory", "<target>", "<baseline>";
                script="benchmark/generators.jl")
=#
using BenchmarkTools

const SUITE = BenchmarkGroup()

SUITE["bimatrix_generators"] = include("bimatrix_generators.jl")
