# GameTheory.jl Benchmarks

This directory contains a benchmark suite in the standard
[BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) format:
[`benchmarks.jl`](benchmarks.jl) defines a `BenchmarkGroup` named `SUITE`,
which can be run standalone or through
[PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl).

Each benchmarked module has its own file, included by `benchmarks.jl`.
Currently covered:

- [`lemke_howson.jl`](lemke_howson.jl): `lemke_howson` and `lemke_howson!`
  (`src/lemke_howson.jl`), under `SUITE["lemke_howson"]`;
- [`support_enumeration.jl`](support_enumeration.jl): `support_enumeration`
  (`src/support_enumeration.jl`), under `SUITE["support_enumeration"]`;
- [`repeated_game.jl`](repeated_game.jl): `outerapproximation`
  (`src/repeated_game.jl`), under `SUITE["repeated_game"]`;
- [`bimatrix_generators.jl`](bimatrix_generators.jl): the game generators
  (`src/generators/bimatrix_generators.jl`), as the separate group
  `GENERATORS_SUITE`, excluded from `SUITE` — run separately; see below.

## What is benchmarked

### `lemke_howson` ([`lemke_howson.jl`](lemke_howson.jl))

The Lemke-Howson algorithm on random `n x n` games (each generated with
its own fixed-seed RNG so that the cases are independent of one another):

| Key | Description |
|:----|:------------|
| `random_n{10,100}` | `lemke_howson` end to end; n = 10 is dominated by setup and allocations, n = 100 by the pivoting steps |
| `random_n{100,200}_capping10` | `lemke_howson` with the capping heuristic of Codenotti et al. (`capping=10`); at n = 200 the plain algorithm can take hundreds of thousands of pivoting steps, which the heuristic avoids |
| `random_n10_prealloc` | `lemke_howson!` with caller-owned output and primary workspace arrays and default keywords (repeated-solve regime) |
| `random_n10_full_workspace` | `lemke_howson!` with the full workspace supplied (allocation floor; no workspace allocations) |

### `support_enumeration` ([`support_enumeration.jl`](support_enumeration.jl))

Support enumeration computes *all* equilibria of a nondegenerate game;
the cost grows combinatorially in the number of actions, so the sizes are
kept small:

| Key | Description |
|:----|:------------|
| `Float/random_n{10,11}` | Random games with `Float64` payoffs (LU-based linear solves) |
| `Rational/identity_n{7,8}` | Identity-matrix games with `Rational{Int}` payoffs (exact-arithmetic path) |

### `repeated_game` ([`repeated_game.jl`](repeated_game.jl))

The outer approximation algorithm of Judd, Yeltekin, and Conklin for the
equilibrium payoff set of a repeated game:

| Key | Description |
|:----|:------------|
| `outerapproximation_nH64` | Prisoner's dilemma with discount factor 0.75, 64 subgradients; the per-iteration linear programs dominate |

### `bimatrix_generators` ([`bimatrix_generators.jl`](bimatrix_generators.jl))

**Run separately**: this subgroup times game construction rather than
equilibrium computation, so it is not part of `SUITE` — whole-suite runs
(standalone or through PkgBenchmark) skip it. The dedicated entry point
[`generators.jl`](generators.jl) exposes it as its `SUITE`; run it
standalone:

```
julia --project=benchmark benchmark/generators.jl
```

or, to run or compare this group with PkgBenchmark, pass the entry point
through the `script` keyword:

```julia
jud = judge("GameTheory", "<target>", "<baseline>";
            script="benchmark/generators.jl")
```

Interactively, `benchmarks.jl` also defines the group as
`GENERATORS_SUITE`; run it the same way as any other subset:

```julia
julia> include("benchmark/benchmarks.jl");

julia> run(GENERATORS_SUITE)
```

Construction of game instances from the test suite of Fearnley, Igwe,
and Savani; the random generators draw a fresh instance per evaluation,
advancing the case's own fixed-seed RNG:

| Key | Description |
|:----|:------------|
| `blotto_game/h{3}_t{62}`, `blotto_game/h{4}_t{21}` | Colonel Blotto games with `(h, t)` hills and troops |
| `ranking_game` | Ranking game with 2000 actions |
| `sgc_game` | SGC game of Sandholm, Gilpin, and Conitzer with k = 500 (4k-1 = 1999 actions per player) |
| `tournament_game` | Tournament game with n = 200, k = 2 |
| `unit_vector_game/avoid_pure_nash_{true,false}` | Unit-vector games with 2000 actions |

## Running the suite standalone

From the repository root, set up the benchmark environment once:

```
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Then run the whole suite (takes a few minutes):

```
julia --project=benchmark benchmark/benchmarks.jl
```

To run interactively, e.g., only a subset:

```julia
julia> include("benchmark/benchmarks.jl");

julia> run(SUITE["lemke_howson"]["random_n100"])
```

## Running with PkgBenchmark.jl

Install PkgBenchmark in your default environment, then, with this package
active (e.g. `julia --project=.`):

```julia
using PkgBenchmark

results = benchmarkpkg("GameTheory")
export_markdown("results.md", results)
```

### Comparing two commits

To evaluate the performance change of a target commit (or branch) relative
to a baseline:

```julia
jud = judge("GameTheory", "<target>", "<baseline>")
export_markdown("judgement.md", jud)
```

For example, to compare the current state of `main` against the previous
commit:

```julia
jud = judge("GameTheory", "main", "main~1")
```

Display the judgment summary:

```julia
julia> show(PkgBenchmark.benchmarkgroup(jud))
```

and the timing estimates of each side:

```julia
julia> show(jud.baseline_results.benchmarkgroup)

julia> show(jud.target_results.benchmarkgroup)
```

Note that `judge` checks out and runs each commit, so uncommitted changes
in the working tree are not included.

For comprehensive usage details, refer to the
[PkgBenchmark documentation](https://juliaci.github.io/PkgBenchmark.jl/stable).
