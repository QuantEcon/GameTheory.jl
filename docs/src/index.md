# GameTheory.jl

[*GameTheory.jl*](https://github.com/QuantEcon/GameTheory.jl) is a [Julia](http://www.julialang.org) package about algorithms and data structures for Game Theory.

## Installation

To install the package, enter the Pkg mode by pressing `]` and run

```julia
add GameTheory
```

## Usage

Once installed, the `GameTheory` package can be used by typing

```@example 1
using GameTheory
```

The Base type `Player` can be created by passing a payoff matrix:

```@example 1
player1 = Player([3 1; 0 2])
```

A 2-player `NormalFormGame` can be created either by passing `Player` instances,

```@example 1
player2 = Player([2 0; 1 3])
g = NormalFormGame((player1, player2))
print(g)
```

or by passing a payoff matrix directly:

```@example 1
payoff_bimatrix = Array{Int}(undef, 2, 2, 2)
payoff_bimatrix[1, 1, :] = [3, 2]
payoff_bimatrix[1, 2, :] = [1, 1]
payoff_bimatrix[2, 1, :] = [0, 0]
payoff_bimatrix[2, 2, :] = [2, 3]
g = NormalFormGame(payoff_bimatrix)
print(g)
```

After constructing a `NormalFormGame`, we can find its Nash Equilibria by using methods of `GameTheory`. For example, `pure_nash` finds all pure action Nash Equilibria by enumeration:

```@example 1
pure_nash(g)
```

Please see the [notebooks](@ref notebooks) on QuantEcon for more details.

## [Notebooks](@id notebooks)

Some notebooks for demonstration are available:

* [Tools for Game Theory](https://nbviewer.jupyter.org/github/QuantEcon/game-theory-notebooks/blob/master/game_theory_jl.ipynb)
* [A Recursive Formulation of Repeated Games](https://nbviewer.jupyter.org/github/QuantEcon/QuantEcon.notebooks/blob/master/recursive_repeated_games.ipynb)

## Library Outline

