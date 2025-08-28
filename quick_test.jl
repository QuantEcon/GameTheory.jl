using Pkg; Pkg.activate(".")
using GameTheory

# Simple test
pd_payoff = [9.0 1.0; 10.0 3.0]
nfg = NormalFormGame(Player(pd_payoff), Player(pd_payoff))
rpd = RepeatedGame(nfg, 0.75)

println("Testing AS algorithm variants...")

# Original method
println("Original:")
r1 = AS(rpd; tol=1e-9, use_optimization=false, incremental_redundancy=false)
println("Size: $(size(r1))")
display(r1)

println("\nOptimized:")
r2 = AS(rpd; tol=1e-9, use_optimization=true, incremental_redundancy=true)
println("Size: $(size(r2))")
display(r2)

# Expected
expected = [3.0 3.0; 3.0 9.75; 9.0 9.0; 9.75 3.0]
println("\nExpected:")
display(expected)