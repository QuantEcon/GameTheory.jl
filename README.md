# GameTheory.jl

[![Build Status](https://github.com/QuantEcon/GameTheory.jl/workflows/CI/badge.svg)](https://github.com/QuantEcon/GameTheory.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/QuantEcon/GameTheory.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/QuantEcon/GameTheory.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://QuantEcon.github.io/GameTheory.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://QuantEcon.github.io/GameTheory.jl/dev)

Algorithms and data structures for game theory in Julia

## Example usage

Create a `NormalFormGame`:

```julia
using GameTheory
player1 = Player([3 3; 2 5; 0 6])
player2 = Player([3 2 3; 2 6 1])
g = NormalFormGame(player1, player2)
println(g)
```
```
3×2 NormalFormGame{2, Int64}:
 [3, 3]  [3, 2]
 [2, 2]  [5, 6]
 [0, 3]  [6, 1]
```

`lrsnash` calls the Nash equilibrium computation routine in [lrslib](http://cgm.cs.mcgill.ca/~avis/C/lrs.html)
(through its Julia wrapper [LRSLib.jl](https://github.com/JuliaPolyhedra/LRSLib.jl)):

```julia
lrsnash(g)
```
```
3-element Vector{Tuple{Vector{Rational{BigInt}}, Vector{Rational{BigInt}}}}:
 ([4//5, 1//5, 0//1], [2//3, 1//3])
 ([0//1, 1//3, 2//3], [1//3, 2//3])
 ([1//1, 0//1, 0//1], [1//1, 0//1])
```

A 2x2x2 `NormalFormGame`:

```julia
g = NormalFormGame((2, 2, 2))
g[1, 1, 1] = [9, 8, 12]
g[2, 2, 1] = [9, 8, 2]
g[1, 2, 2] = [3, 4, 6]
g[2, 1, 2] = [3, 4, 4]
println(g)
```
```
2×2×2 NormalFormGame{3, Float64}:
[:, :, 1] =
 [9.0, 8.0, 12.0]  [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]   [9.0, 8.0, 2.0]

[:, :, 2] =
 [0.0, 0.0, 0.0]  [3.0, 4.0, 6.0]
 [3.0, 4.0, 4.0]  [0.0, 0.0, 0.0]
```

`hc_solve` computes all isolated Nash equilibria of an N-player game by using
[HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl):

```julia
NEs = hc_solve(g)
```
```
9-element Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}:
 ([2.63311e-36, 1.0], [0.333333, 0.666667], [0.333333, 0.666667])
 ([0.25, 0.75], [1.0, 0.0], [0.25, 0.75])
 ([0.0, 1.0], [0.0, 1.0], [1.0, 0.0])
 ([0.25, 0.75], [0.5, 0.5], [0.333333, 0.666667])
 ([0.5, 0.5], [0.5, 0.5], [1.0, 1.37753e-40])
 ([1.0, 0.0], [0.0, 1.0], [0.0, 1.0])
 ([0.5, 0.5], [0.333333, 0.666667], [0.25, 0.75])
 ([1.0, 0.0], [1.0, 9.40395e-38], [1.0, -9.40395e-38])
 ([0.0, 1.0], [1.0, 0.0], [0.0, 1.0])
```

See the tutorials for further examples.

## Implemented algorithms

### Nash equilibrium computation

* [`pure_nash`](https://quantecon.github.io/GameTheory.jl/stable/lib/computing_nash_equilibria.html#GameTheory.pure_nash-Tuple{NormalFormGame}):
  Find all pure-action Nash equilibria of an N-player game (if any)
* [`support_enumeration`](https://quantecon.github.io/GameTheory.jl/stable/lib/computing_nash_equilibria.html#GameTheory.support_enumeration-Union{Tuple{NormalFormGame{2,%20T}},%20Tuple{T}}%20where%20T):
  Find all mixed-action Nash equilibria of a two-player nondegenerate game
* [`lrsnash`](https://quantecon.github.io/GameTheory.jl/stable/lib/computing_nash_equilibria.html#GameTheory.lrsnash-Tuple{NormalFormGame{2,%20%3C:Union{Int64,%20Rational}}}):
  Find all mixed-action Nash equilibria (or equilibrium components) of a two-player game
* [`hc_solve`](https://quantecon.github.io/GameTheory.jl/stable/lib/computing_nash_equilibria.html#GameTheory.hc_solve-Union{Tuple{NormalFormGame{N}},%20Tuple{N}}%20where%20N):
  Find all isolated mixed-action Nash equilibria of an N-player game

### Learning/evolutionary dynamics

* [`BRD`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.BRD):
  Best response dynamics
* [`KMR`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.KMR):
  Best response with mutations dynamics of Kandori-Mailath-Rob
* [`SamplingBRD`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.SamplingBRD):
  Sampling best response dynamics
* [`FictitiousPlay`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.FictitiousPlay):
  Fictitious play
* [`StochasticFictitiousPlay`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.StochasticFictitiousPlay):
  Stochastic fictitious play
* [`LocalInteraction`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.LocalInteraction):
  Local interaction dynamics
* [`LogitDynamics`](https://quantecon.github.io/GameTheory.jl/stable/lib/learning_algorithms.html#GameTheory.LogitDynamics):
  Logit dynamics

### Repeated games

* [`outerapproximation`](https://quantecon.github.io/GameTheory.jl/stable/lib/repeated_games.html#GameTheory.outerapproximation-Tuple{RepeatedGame{2}}):
  Equilibrium payoff computation algorithm by Judd-Yeltekin-Conklin
* [`AS`](https://quantecon.github.io/GameTheory.jl/stable/lib/repeated_games.html#GameTheory.AS-Union{Tuple{RepeatedGame{2,%20T,%20TD}},%20Tuple{TD},%20Tuple{T}}%20where%20{T,%20TD}):
  Equilibrium payoff computation algorithm by Abreu-Sannikov

## Tutorials

* [Tools for Game Theory in GameTheory.jl](https://nbviewer.org/github/QuantEcon/game-theory-notebooks/blob/main/game_theory_jl.ipynb)
* [A Recursive Formulation of Repeated Games](https://nbviewer.org/github/QuantEcon/QuantEcon.notebooks/blob/master/recursive_repeated_games.ipynb)
* [Abreu-Sannikov Algorithm for Repeated Two-player Games](https://nbviewer.jupyter.org/github/QuantEcon/game-theory-notebooks/blob/main/as_jl.ipynb)

See also the [`game_theory`](https://quanteconpy.readthedocs.io/en/latest/game_theory.html) submodule of
[`QuantEcon.py`](https://github.com/QuantEcon/QuantEcon.py),
the Python counterpart of this package.
