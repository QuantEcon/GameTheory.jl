# GameTheory.jl Benchmarks

## Running the benchmark suite

Use [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl):

```jl
using PkgBenchmark
```

As an example, let us run the benchmarks (defined in `benchmarks.jl`) and compare the performance changes for the two commits
[`32e6090`](https://github.com/QuantEcon/GameTheory.jl/commit/32e60906bdf34f39bc535fc1235e5da5e261d1c4) (target)
and
[`b5031c3`](https://github.com/QuantEcon/GameTheory.jl/commit/b5031c3cf24e41b124e91f4881c9a543eed19ecc) (baseline):

```jl
jud = judge("GameTheory", "32e6090", "b5031c3")
```

To show the results:

```jl
julia> show(PkgBenchmark.benchmarkgroup(jud))
2-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
    tags: ["support_enumeration"]
    "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (8, 8) => TrialJudgement(+5.42% => regression)
      (7, 7) => TrialJudgement(+5.27% => regression)
    "Float" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (10, 10) => TrialJudgement(-1.00% => invariant)
      (11, 11) => TrialJudgement(-0.25% => invariant)
...
```

To show the timing estimates for the baseline:

```jl
julia> show(jud.baseline_results.benchmarkgroup)
2-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
    tags: ["support_enumeration"]
    "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (8, 8) => Trial(128.530 ms)
      (7, 7) => Trial(29.616 ms)
    "Float" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (10, 10) => Trial(96.379 ms)
      (11, 11) => Trial(410.365 ms)
...
```

and for the target:

```jl
julia> show(jud.target_results.benchmarkgroup)
2-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
    tags: ["support_enumeration"]
    "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (8, 8) => Trial(135.496 ms)
      (7, 7) => Trial(31.177 ms)
    "Float" => 2-element BenchmarkTools.BenchmarkGroup:
      tags: []
      (10, 10) => Trial(95.419 ms)
      (11, 11) => Trial(409.337 ms)
...
```

To run a script file other than `benchmarks.jl`:

```jl
using GameTheory
results = benchmarkpkg(
    "GameTheory",
    script="$(dirname(pathof(GameTheory)))/../benchmark/generators.jl"
)
```

For more usage information, see the [PkgBenchmark documentation](https://juliaci.github.io/PkgBenchmark.jl/stable).
