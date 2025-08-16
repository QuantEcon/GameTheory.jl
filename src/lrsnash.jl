using LRSLib: buildrep, solve_nash

"""
    lrsnash(g)

Compute in exact arithmetic all extreme mixed-action Nash equilibria of
a 2-player normal form game with Integer or Rational payoffs. This
function calls the Nash equilibrium computation routine in `lrslib`
(through its Julia wrapper `LRSLib.jl`) which is based on the
"lexicographic reverse search" vertex enumeration algorithm.

# Arguments

- `g::NormalFormGame{2,<:RatOrInt}`: 2-player NormalFormGame instance
  with Integer or Rational payoffs.

# Returns

- `::Vector{NTuple{2,Vector{Rational{BigInt}}}}`: Vector of mixed-action
  Nash equilibria.

# Examples

A degenerate game example:

```julia
julia> player1 = Player([3 3; 2 5; 0 6]);

julia> player2 = Player([3 2 3; 3 6 1]);

julia> g = NormalFormGame(player1, player2);

julia> println(g)
3×2 NormalFormGame{2, Int64}:
 [3, 3]  [3, 3]
 [2, 2]  [5, 6]
 [0, 3]  [6, 1]

julia> lrsnash(g)
3-element Vector{Tuple{Vector{Rational{BigInt}}, Vector{Rational{BigInt}}}}:
 ([1//1, 0//1, 0//1], [1//1, 0//1])
 ([1//1, 0//1, 0//1], [2//3, 1//3])
 ([0//1, 1//3, 2//3], [1//3, 2//3])
```

The set of Nash equilibria of this degenerate game consists of an
isolated equilibrium, the third output, and a non-singleton equilibrium
component, the extreme points of which are given by the first
two outputs.

# References

- D. Avis, G. Rosenberg, R. Savani, and B. von Stengel, "Enumeration of
  Nash Equilibria for Two-Player Games," Economic Theory (2010), 9-37.
"""
function lrsnash(g::NormalFormGame{2,<:RatOrInt})
    hrs = [buildrep(i, g.players[3-i].payoff_array) for i in 1:2]
    NEs = solve_nash(hrs...)
    return NEs
end
