# Copilot Instructions for GameTheory.jl

## Repository Overview

GameTheory.jl is a Julia package that implements algorithms and data structures for game theory. The package provides tools for:

- **Normal Form Games**: Creating and analyzing strategic games with multiple players
- **Nash Equilibrium Computation**: Finding pure and mixed strategy Nash equilibria
- **Learning/Evolutionary Dynamics**: Simulating how players' strategies evolve over time
- **Repeated Games**: Analyzing games played repeatedly over time

## Working Effectively

### Setup and tests:
- Install dependencies: `julia --project=. -e "using Pkg; Pkg.instantiate()"`
- Run all tests: `julia --project=. -e "using Pkg; Pkg.test()"`
- Run a single test file: `julia --project=. -e 'using GameTheory, Test; include("test/util.jl"); include("test/test_pure_nash.jl")'`. The test files are written to be included from `test/runtests.jl`, which always includes all of them: a test file cannot be run on its own (`julia test/test_pure_nash.jl` fails), since `GameTheory`, `Test` and the helpers in `test/util.jl` are loaded there. `test/test_aqua.jl` needs Aqua.jl, a test-only dependency, so it runs only through `Pkg.test()`.
- Julia 1.10+ is required (the `[compat]` bound in `Project.toml` is `julia = "1.10"`).

### Documentation:
- Build with the local checkout (see `docs/README.md`): `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`
- CI also runs the doctests in the docstrings: `julia --project=docs -e 'using Documenter: doctest; using GameTheory; doctest(GameTheory)'`. Run them after changing a docstring example.

### Benchmarks:
- The benchmark suite under `benchmark/` is in the standard BenchmarkTools.jl format (see `benchmark/README.md` for full usage). Setup (once): `julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`; run: `julia --project=benchmark benchmark/benchmarks.jl`.
- Compare two commits with PkgBenchmark.jl: `judge("GameTheory", "<target>", "<baseline>")`; `judge` runs committed states, so commit your changes first.
- Benchmarks are not run in CI. For performance-sensitive changes, report before/after numbers from this suite in the PR description.

## Key Concepts and Types

### Core Types
- `Player{N,T}`: Represents a player with an N-dimensional payoff array of type T
- `NormalFormGame{N,T}`: Represents an N-player normal form game
- `RepeatedGame`: For repeated games analysis

### Type Aliases
- `PureAction = Integer`: A pure strategy (single action choice)
- `MixedAction{T} = Vector{T}`: A mixed strategy (probability distribution over actions)
- `Action{T} = Union{PureAction,MixedAction{T}}`: Either pure or mixed action
- `ActionProfile{N,T}`: A tuple of N actions (one per player)

### Important Constants
- `RatOrInt = Union{Rational,Integer}`: Used for exact arithmetic in Nash equilibrium computations

## Code Organization

### Core Modules (`src/`)
- `normal_form_game.jl`: Main game representation and basic operations
- `pure_nash.jl`: Pure strategy Nash equilibrium computation
- `support_enumeration.jl`: Mixed strategy Nash equilibria via support enumeration
- `vertex_enumeration.jl`: Mixed strategy Nash equilibria via vertex enumeration
- `lemke_howson.jl`: A mixed strategy Nash equilibrium by the Lemke-Howson algorithm
- `lrsnash.jl`: Nash equilibria using LRS library (vertex enumeration)
- `homotopy_continuation.jl`: Nash equilibria using polynomial homotopy continuation
- `repeated_game.jl`: Tools for repeated games analysis
- `random.jl`: Random game generation utilities
- `game_converters.jl`: Readers and writers for the GameTracer `.gam` and Gambit `.nfg` formats
- `util.jl`: General utility functions

### Learning Algorithms (`src/`)
- `fictplay.jl`: Fictitious play dynamics
- `localint.jl`: Local interaction dynamics on networks
- `brd.jl`: Best response dynamics and variants (BRD, KMR, SamplingBRD)
- `logitdyn.jl`: Logit choice dynamics

### Generators (`src/generators/`)
- `Generators.jl`: The `Generators` submodule
- `bimatrix_generators.jl`: Generators of 2-player games

## Development Guidelines

### Workspace and Allocation Conventions
- Mutating solver variants (e.g. `lemke_howson!`) follow the QuantEcon.jl `lcp_lemke!` argument layout: caller-owned output and primary workspace arrays are positional, auxiliary workspace arrays are keyword arguments defaulting to `nothing` and materialized lazily (keyword defaults are evaluated at call time even if the function returns early)
- When validating caller-supplied arguments in exported functions, prefer throwing exceptions to `@assert`: `DimensionMismatch` for array-size mismatches, `ArgumentError` for invalid values or aliasing
- Check non-aliasing with `Base.mightalias` among caller-supplied arrays that the algorithm overwrites before reading others (e.g. `argmins` vs `bases` in `lemke_howson!`), and between the members of same-shaped array pairs
- In allocation tests, do not assert `@allocated ... == 0` or hard-code byte counts: older Julia versions heap-allocate returned non-isbits immutable objects that newer versions elide; bound the measurement by a baseline measured in the same escape pattern (see `test/test_lemke_howson.jl`), and scope "allocation-free" claims to machine-float element types

### Testing
- Most source files have a corresponding test file (e.g., `src/normal_form_game.jl` → `test/test_normal_form_game.jl`), included from `test/runtests.jl`; the tests for `src/generators/` are in `test/generators/` and are included from `test/generators/runtests.jl`. Add new test files to the respective `runtests.jl`.
- Each test file loads the packages it uses, other than `GameTheory` and `Test`, with its own `using` statements, so that it does not depend on the files included before it.

## Common Patterns and Idioms

### Game Creation
```julia
# Create players with payoff matrices
player1 = Player([3 0; 5 1])  # 2x2 payoff matrix
player2 = Player([3 5; 0 1])
game = NormalFormGame(player1, player2)

# Or create directly from an array of payoff tuples, one tuple per action
# profile: payoffs[a1, a2] = (player 1's payoff, player 2's payoff)
game = NormalFormGame([(3, 3) (0, 0); (5, 5) (1, 1)])

# Note: there is no NormalFormGame(matrix1, matrix2) constructor; a single
# square matrix NormalFormGame(A) constructs a *symmetric* 2-player game
```

### Action Handling
```julia
# Pure actions are integers (1-indexed)
pure_action = 1

# Mixed actions are probability vectors
mixed_action = [0.6, 0.4]  # 60% action 1, 40% action 2

# Action profiles for multiple players
profile = (1, 2)  # Player 1 plays action 1, Player 2 plays action 2
```

### Nash Equilibrium Computation
```julia
# Different methods for different game types
pure_equilibria = pure_nash(game)
mixed_equilibria = support_enumeration(game)  # 2-player only
all_equilibria = lrsnash(game)  # Exact rational arithmetic
```

### Learning Dynamics
```julia
# Set up dynamics
dynamics = FictitiousPlay(game)
initial_actions = (1, 1)

# Simulate
final_actions = play(dynamics, initial_actions, num_reps=1000)
history = time_series(dynamics, 100, initial_actions)
```

## Contribution Conventions

### Keeping these instructions up to date:
- This file (`.github/copilot-instructions.md`) is the single source of repository instructions for AI agents (`AGENTS.md` only points here). Update it whenever necessary as part of the change that makes it outdated: when adding or restructuring files, changing workflows or conventions, bumping the required Julia version, or learning a repository-specific pitfall worth passing on. Stale instructions are worse than none.

### Commit, PR and issue titles:
- Prefix commit subjects and PR and issue titles with the change type: `ENH:` for new features and enhancements, `FIX:` for bug fixes (use `FIX:`, not `BUG:`), `PERF:` for performance changes, `TEST:` for test-only changes, `DOC:` for documentation, `MAINT:` for maintenance, `RFC:` for refactoring, `CI:` for CI configuration.

### PR descriptions and other GitHub text:
- Do not hard-wrap lines: keep each paragraph and each bullet point on a single line (GitHub soft-wraps; manual line breaks render poorly).
- Wrap every `@`-prefixed token in backticks, in PR and issue descriptions, comments, and commit messages alike: a bare `@name` (typically a Julia macro such as `@inferred` or `@test`) is a GitHub mention and notifies whoever owns that username.
- State explicitly whether the change is behavior-preserving; if it removes or changes anything that appears in the published API documentation (including docstrings picked up by `@autodocs`), flag it as a breaking change in the PR description. (Release notes are created on the [GitHub Releases page](https://github.com/QuantEcon/GameTheory.jl/releases) at release time.)
- For performance PRs, include measured before/after numbers from the benchmark suite.

### Declaration of AI assistance:
- End every PR description, issue, and comment written with AI assistance with a one-line declaration naming the tool and the model, in the form `Assisted-by: <tool> (<model>)` or `Generated with <tool> (<model>)`; for example `Assisted-by: Claude Code (Claude Fable 5.1)` or `🤖 Generated with [Claude Code](https://claude.com/claude-code) (Claude Fable 5.1)`.
- In commit messages, keep the `Co-Authored-By:` trailer that the tool emits.

### Cross-language parity with QuantEcon.py:
- This package is the Julia counterpart of the [`game_theory`](https://quanteconpy.readthedocs.io/en/latest/game_theory.html) submodule of [QuantEcon.py](https://github.com/QuantEcon/QuantEcon.py); for example, `src/vertex_enumeration.jl` mirrors `quantecon/game_theory/vertex_enumeration.py`.
- A bug found in one implementation likely exists in the other: check the sibling and fix (or at least report) both.
- Behavioral changes should keep the two implementations consistent; note the corresponding PR/issue of the sibling repository in the PR description.

## Validation

After making code changes, run the test files relevant to the change, then the full test suite (`Pkg.test()`) before finalizing. If a test fails, investigate and fix before committing.

## CI

- `.github/workflows/ci.yml` runs the tests on the latest stable Julia 1.x on Ubuntu, Windows, and macOS, runs the doctests, and builds the documentation.
- `.github/workflows/ci-nighty.yml` runs the tests on Julia nightly on Ubuntu.

## Resources

- [Package Documentation](https://quantecon.github.io/GameTheory.jl/stable/)
- [QuantEcon Lectures](https://julia.quantecon.org/) for economic applications
- [Python version (QuantEcon.py game_theory submodule)](https://quanteconpy.readthedocs.io/en/latest/game_theory.html)