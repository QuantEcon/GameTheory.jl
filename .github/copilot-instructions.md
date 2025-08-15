# Copilot Instructions for GameTheory.jl

## Repository Overview

GameTheory.jl is a Julia package that implements algorithms and data structures for game theory. The package provides tools for:

- **Normal Form Games**: Creating and analyzing strategic games with multiple players
- **Nash Equilibrium Computation**: Finding pure and mixed strategy Nash equilibria
- **Learning/Evolutionary Dynamics**: Simulating how players' strategies evolve over time
- **Repeated Games**: Analyzing games played repeatedly over time

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
- `RatOrInt = Union{Rational,Int}`: Used for exact arithmetic in Nash equilibrium computations

## Code Organization

### Core Modules (`src/`)
- `normal_form_game.jl`: Main game representation and basic operations
- `pure_nash.jl`: Pure strategy Nash equilibrium computation
- `support_enumeration.jl`: Mixed strategy Nash equilibria via support enumeration
- `lrsnash.jl`: Nash equilibria using LRS library (vertex enumeration)
- `homotopy_continuation.jl`: Nash equilibria using polynomial homotopy continuation
- `repeated_game.jl`: Tools for repeated games analysis
- `random.jl`: Random game generation utilities
- `util.jl`: General utility functions

### Learning Algorithms (`src/`)
- `fictplay.jl`: Fictitious play dynamics
- `localint.jl`: Local interaction dynamics on networks
- `brd.jl`: Best response dynamics and variants (BRD, KMR, SamplingBRD)
- `logitdyn.jl`: Logit choice dynamics

### Generators (`src/generators/`)
- `Generators.jl`: Game generation utilities
- Various specialized game generators (bimatrix, coordination, etc.)

## Development Guidelines

### Julia Conventions
- Follow standard Julia naming conventions (lowercase with underscores for functions, CamelCase for types)
- Use `!` suffix for mutating functions (e.g., `play!`)
- Provide both mutating and non-mutating versions where appropriate
- Use multiple dispatch effectively for different input types

### Mathematical Precision
- Use `Rational` types for exact computations when needed (especially Nash equilibria)
- Be careful with floating-point comparisons; use appropriate tolerances
- Prefer numerically stable algorithms for matrix operations

### Performance Considerations
- Use type-stable functions with proper type annotations
- Leverage Julia's multiple dispatch for efficient specialization
- Consider using `@inbounds` for performance-critical loops (with proper bounds checking)
- Use views (`@view`) instead of creating temporary arrays when possible

### Testing Patterns
- Each source file has a corresponding test file (e.g., `src/normal_form_game.jl` → `test/test_normal_form_game.jl`)
- Use `@testset` to group related tests
- Include edge cases and error conditions in tests
- Test both pure and mixed strategy scenarios
- Verify mathematical properties (e.g., Nash equilibrium conditions)

### Documentation Standards
- Use comprehensive docstrings with examples
- Include mathematical background for algorithms
- Document function arguments with types and descriptions
- Provide usage examples for complex functions
- Reference academic papers for algorithm implementations

## Key Dependencies and Their Uses

### Mathematical Libraries
- `LinearAlgebra`: Matrix operations for payoff computations
- `Distributions`: Probability distributions for stochastic processes
- `Combinatorics`: Generating combinations/permutations of strategies
- `StaticArrays`: High-performance small arrays (where applicable)

### Optimization Libraries
- `MathOptInterface` (MOI): Abstract interface for optimization problems
- `HiGHS`: High-performance linear/mixed-integer programming solver
- `Polyhedra`: Geometric computations with polytopes
- `CDDLib`: Convex hull and vertex enumeration
- `LRSLib`: Linear reverse search for vertex/facet enumeration

### Specialized Libraries
- `HomotopyContinuation`: Polynomial system solving for Nash equilibria
- `QuantEcon`: Economic and quantitative modeling utilities
- `Graphs`: Network analysis for local interaction models

## Common Patterns and Idioms

### Game Creation
```julia
# Create players with payoff matrices
player1 = Player([3 0; 5 1])  # 2x2 payoff matrix
player2 = Player([3 5; 0 1])
game = NormalFormGame(player1, player2)

# Or create directly
game = NormalFormGame([3 0; 5 1], [3 5; 0 1])
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

## Testing and Quality Assurance

### Running Tests
- Use `julia --project=. -e "using Pkg; Pkg.test()"` to run full test suite
- Individual test files can be run with `julia test/test_filename.jl`
- CI runs tests on Julia 1.x across Linux, Windows, and macOS

### Code Quality
- Follow the existing code style in the repository
- Add tests for new functionality
- Update documentation for user-facing changes
- Consider numerical stability and edge cases
- Test with both small and large games where applicable

## Mathematical Background

When implementing new algorithms, ensure understanding of:
- **Game Theory**: Nash equilibria, dominant strategies, Pareto efficiency
- **Linear Algebra**: Matrix operations, eigenvalues, linear systems
- **Optimization**: Linear programming, convex optimization
- **Probability**: Mixed strategies as probability distributions
- **Numerical Methods**: Stability, convergence, precision issues

## Common Gotchas

1. **Indexing**: Julia uses 1-based indexing, not 0-based
2. **Mutability**: Be careful with shared references to arrays
3. **Type Stability**: Avoid changing variable types within functions
4. **Broadcasting**: Use `.` for element-wise operations (`.*`, `.+`, etc.)
5. **Rational Arithmetic**: When using `Rational` types, operations may be slower but exact
6. **Memory Management**: Large games can consume significant memory; consider algorithmic complexity

## Resources

- [Julia Documentation](https://docs.julialang.org/)
- [Game Theory Textbooks]: Fudenberg & Tirole, Myerson, etc.
- [Package Documentation](https://quantecon.github.io/GameTheory.jl/stable/)
- [QuantEcon Lectures](https://lectures.quantecon.org/) for economic applications
- [Python version (QuantEcon.py game_theory submodule)](https://quanteconpy.readthedocs.io/en/latest/game_theory.html)