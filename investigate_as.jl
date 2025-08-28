#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

using GameTheory
using Printf

function investigate_as_results()
    println("🔍 AS Algorithm Result Investigation")
    println("=" ^ 60)
    
    # Create test case
    pd_payoff = [9.0 1.0; 10.0 3.0]
    nfg = NormalFormGame(Player(pd_payoff), Player(pd_payoff))
    rpd = RepeatedGame(nfg, 0.75)
    
    println("Test Game: 2x2 Prisoner's Dilemma, δ=0.75")
    println()
    
    # Get results from both methods
    println("Getting results from original method...")
    result1 = AS(rpd; tol=1e-9, use_optimization=false, incremental_redundancy=false, verbose=true)
    println("Original result:")
    display(result1)
    println()
    
    println("Getting results from optimized method...")
    result2 = AS(rpd; tol=1e-9, use_optimization=true, incremental_redundancy=true, verbose=true)
    println("Optimized result:")
    display(result2)
    println()
    
    # Expected vertices
    expected = [3.0 3.0; 3.0 9.75; 9.0 9.0; 9.75 3.0]
    println("Expected vertices:")
    display(expected)
    println()
    
    function detailed_check(vertices, name)
        println("$name detailed check:")
        for i in 1:size(vertices, 1)
            println(@sprintf("  Vertex %d: [%.6f, %.6f]", i, vertices[i,1], vertices[i,2]))
        end
        
        println("  Matching against expected:")
        matches = 0
        for i in 1:size(expected, 1)
            found = false
            for j in 1:size(vertices, 1)
                if abs(vertices[j,1] - expected[i,1]) < 1e-6 && abs(vertices[j,2] - expected[i,2]) < 1e-6
                    found = true
                    matches += 1
                    println(@sprintf("    ✅ Expected [%.2f, %.2f] matches vertex %d", expected[i,1], expected[i,2], j))
                    break
                end
            end
            if !found
                println(@sprintf("    ❌ Expected [%.2f, %.2f] NOT found", expected[i,1], expected[i,2]))
            end
        end
        println(@sprintf("    Total matches: %d/%d", matches, size(expected, 1)))
        println()
        return matches == size(expected, 1)
    end
    
    correct1 = detailed_check(result1, "Original")
    correct2 = detailed_check(result2, "Optimized")
    
    if !correct2
        println("⚠️  Correctness issue with optimized method detected!")
        println("This suggests the optimizations may have bugs.")
        
        # Try with just one optimization at a time
        println("\nTesting optimizations individually:")
        
        println("Testing use_optimization=true, incremental_redundancy=false...")
        result3 = AS(rpd; tol=1e-9, use_optimization=true, incremental_redundancy=false)
        println("Result:")
        display(result3)
        correct3 = detailed_check(result3, "Optimization only")
        
        println("Testing use_optimization=false, incremental_redundancy=true...")
        result4 = AS(rpd; tol=1e-9, use_optimization=false, incremental_redundancy=true)
        println("Result:")
        display(result4)
        correct4 = detailed_check(result4, "Incremental only")
        
        if correct3 && !correct4
            println("🔍 Issue appears to be with incremental_redundancy=true")
        elif !correct3 && correct4
            println("🔍 Issue appears to be with use_optimization=true")
        elif !correct3 && !correct4
            println("🔍 Issues with both optimizations")
        else
            println("🤔 Issue only appears when both are combined")
        end
    else
        println("✅ Both methods produce correct results!")
    end
end

try
    investigate_as_results()
catch e
    println("Error: $e")
    Base.show_backtrace(stdout, catch_backtrace())
end