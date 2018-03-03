# Games.jl Benchmarks

## Running the benchmark suite

Use [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl):

```jl
using PkgBenchmark
```

As an example, let us run the benchmarks (defined in `benchmarks.jl`) and compare the performance changes for the two commits
[`68f3c4b`](https://github.com/QuantEcon/Games.jl/commit/68f3c4bef03554a00384350a047f1e95abd865df) (target)
and
[`d6682de`](https://github.com/QuantEcon/Games.jl/commit/d6682deb9fdae6f16a89b17fbeee9061d763710a) (baseline):

```jl
jud = judge("Games", "68f3c4b", "d6682de")
```

To show the results:

```jl
julia> showall(PkgBenchmark.benchmarkgroup(jud))
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
	  tags: ["support_enumeration"]
	  "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(8, 8)" => TrialJudgement(-1.86% => invariant)
		  "(7, 7)" => TrialJudgement(-1.68% => invariant)
	  "Float" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(10, 10)" => TrialJudgement(-8.62% => improvement)
		  "(11, 11)" => TrialJudgement(-10.96% => improvement)
```

To show the timing estimates for the baseline:

```jl
julia> showall(jud.baseline_results.benchmarkgroup)
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
	  tags: ["support_enumeration"]
	  "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(8, 8)" => Trial(442.854 ms)
		  "(7, 7)" => Trial(97.427 ms)
	  "Float" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(10, 10)" => Trial(207.076 ms)
		  "(11, 11)" => Trial(867.757 ms)
```

and for the target:

```jl
julia> showall(jud.target_results.benchmarkgroup)
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "support_enumeration" => 2-element BenchmarkTools.BenchmarkGroup:
	  tags: ["support_enumeration"]
	  "Rational" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(8, 8)" => Trial(434.614 ms)
		  "(7, 7)" => Trial(95.789 ms)
	  "Float" => 2-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "(10, 10)" => Trial(189.223 ms)
		  "(11, 11)" => Trial(772.689 ms)
```

For more usage information, see the [PkgBenchmark documentation](https://juliaci.github.io/PkgBenchmark.jl/stable).
