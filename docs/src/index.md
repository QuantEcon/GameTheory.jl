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

### Creating a game

The Base type `Player` can be created by passing a payoff matrix:

```@example 1
player1 = Player([3 1; 0 2])
```

Here the rows of the payoff matrix correspond to player 1's own actions and
the columns to the opponent's actions.

A 2-player `NormalFormGame` can be created either by passing `Player` instances,

```@example 1
player2 = Player([2 0; 1 3])
g = NormalFormGame((player1, player2))
```

or by passing an array of tuples representing payoff profiles:

```@example 1
g = NormalFormGame([(3, 2) (1, 1)
                    (0, 0) (2, 3)])
```

or by passing a payoff matrix directly:

```@example 1
payoff_bimatrix = Array{Int}(undef, 2, 2, 2)
payoff_bimatrix[1, 1, :] = [3, 2]
payoff_bimatrix[1, 2, :] = [1, 1]
payoff_bimatrix[2, 1, :] = [0, 0]
payoff_bimatrix[2, 2, :] = [2, 3]
g = NormalFormGame(payoff_bimatrix)
```

### Payoff array conventions

Each player's `payoff_array` is indexed with the player's *own action first*:
for player `i` in an N-player game, the first axis corresponds to player
`i`'s own actions, and the `j`-th axis, `j = 2, ..., N`, to the actions of
player `i+j-1` (modulo `N`). In the 2-player game `g` constructed above,

```@example 1
player1.payoff_array[1, 2] == 1 && player2.payoff_array[2, 1] == 1
```

both give the payoffs under the action profile in which player 1 plays
action 1 and player 2 plays action 2 — each player's array is indexed with
that player's own action first. Similarly, in a 3-player game,
`g.players[2].payoff_array[a2, a3, a1]` is player 2's payoff under the
action profile `(a1, a2, a3)`.

### Computing Nash equilibria

After constructing a `NormalFormGame`, we can find its Nash equilibria by
using methods of `GameTheory`. For example, `pure_nash` finds all pure-action
Nash equilibria by enumeration:

```@example 1
pure_nash(g)
```

The game also has a mixed-action Nash equilibrium: `vertex_enumeration`
finds all Nash equilibria of a two-player nondegenerate game, pure and mixed:

```@example 1
vertex_enumeration(g)
```

See [Computing Nash Equilibria](@ref computing_nash_equilibria) for the other
solvers available.

## [Notebooks](@id notebooks)

Some notebooks for demonstration are available:

* [Tools for Game Theory](https://nbviewer.jupyter.org/github/QuantEcon/game-theory-notebooks/blob/main/game_theory_jl.ipynb)
* [A Recursive Formulation of Repeated Games](https://nbviewer.jupyter.org/github/QuantEcon/QuantEcon.notebooks/blob/main/recursive_repeated_games.ipynb)

## Library Outline

* [Base Types and Methods](@ref base_types_and_methods)

* [Game Generators](@ref game_generators)

* [Computing Nash Equilibria](@ref computing_nash_equilibria)

* [Learning Algorithms](@ref learning_algorithms)

* [Repeated Games](@ref repeated_games)

* [Utilities](@ref util)

