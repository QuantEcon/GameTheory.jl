# Games.jl

*Games.jl* is a [`Julia`](http://www.julialang.org) package for computation of Game Theory.

## Installation

*Games.jl* is an unregistered package that currently under development.

To install the package, open a `Julia` session and type

```julia
Pkg.clone("https://github.com/QuantEcon/Games.jl")
```

## Usage

Once installed, the `Games` package can be used by typing

```@example 1
using Games
```

The Base type `Player` can be created by passing a payoff matrix.

```@example 1
player1 = Player([3 1; 0 2])
```

A 2-player `NormalFormGame` can be created either by passing `Player` instances.

```@example 1
player2 = Player([2 0; 1 3])
g = NormalFormGame((player1, player2))
```

or by passing a payoff bimatrix directly.

```@example 1
payoff_bimatrix = Array(Int, 2, 2, 2)
payoff_bimatrix[1, 1, :] = [3, 2]
payoff_bimatrix[1, 2, :] = [1, 1]
payoff_bimatrix[2, 1, :] = [0, 0]
payoff_bimatrix[2, 2, :] = [2, 3]
g = NormalFormGame(payoff_bimatrix)
```

After construction of `NormalFormGame`, we can find its Nash Equilibria by using methods of `Games`, for example, `pure_nash` finds all pure action Nash Equilibria by enumeration.

```@example 1
pure_nash(g)
```

Please see [lectures](@ref lectures) on QuantEcon for more details.

## [Lectures](@id lectures)

Some lectures of this package is available on [QuantEcon](https://lectures.quantecon.org):

* [Tools for Game Theory](http://nbviewer.jupyter.org/github/QuantEcon/QuantEcon.notebooks/blob/master/game_theory_jl.ipynb)

* [A Recursive Formulation of Repeated Games](https://nbviewer.jupyter.org/github/QuantEcon/QuantEcon.notebooks/blob/master/recursive_repeated_games.ipynb)

## Library Outline

