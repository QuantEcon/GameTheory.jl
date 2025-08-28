#!/usr/bin/env julia

"""
Minimal Performance Test for AS Algorithm Optimizations

This directly tests the AS function with and without optimizations
to empirically measure the performance difference between versions.
"""

using Pkg
Pkg.activate(".")

# Load the package
using GameTheory
using Printf

function simple_time(f, n=3)
    """Simple timing function - run n times, return minimum"""
    times = Float64[]
    for i in 1:n
        GC.gc()  # Force garbage collection
        t = @elapsed f()
        push!(times, t)
    end
    return minimum(times)
end

function run_benchmark()
    println("🔬 AS Algorithm Performance Test")
    println("=" ^ 50)
    println("Comparing optimization flags from PR #65:")
    println("- Original: use_optimization=false, incremental_redundancy=false")
    println("- Optimized: use_optimization=true, incremental_redundancy=true")
    println()

    # Create a test case
    pd_payoff = [9.0 1.0; 10.0 3.0]
    nfg = NormalFormGame(Player(pd_payoff), Player(pd_payoff))
    rpd = RepeatedGame(nfg, 0.75)
    
    println("Test Game: 2x2 Prisoner's Dilemma, δ=0.75")
    println()
    
    # Test original method
    println("Testing Original Method...")
    original_time = simple_time(() -> AS(rpd; tol=1e-9, 
                                        use_optimization=false, 
                                        incremental_redundancy=false))
    println(@sprintf("  Time: %.4f seconds", original_time))
    
    # Test optimized method  
    println("Testing Optimized Method...")
    optimized_time = simple_time(() -> AS(rpd; tol=1e-9,
                                         use_optimization=true, 
                                         incremental_redundancy=true))
    println(@sprintf("  Time: %.4f seconds", optimized_time))
    
    # Calculate speedup
    speedup = original_time / optimized_time
    
    println()
    println("=" ^ 50)
    println("RESULTS:")
    println(@sprintf("Original time:   %.4f s", original_time))
    println(@sprintf("Optimized time:  %.4f s", optimized_time))
    println(@sprintf("Speedup:         %.2fx", speedup))
    
    if speedup > 1.1
        println("✅ Significant improvement - optimizations are effective!")
    elseif speedup > 1.0
        println("✅ Modest improvement - optimizations help")
    elseif speedup > 0.9
        println("➡️  Similar performance")
    else
        println("📉 Optimizations appear slower")
    end
    
    # Verify correctness
    println()
    println("Checking correctness...")
    result1 = AS(rpd; tol=1e-9, use_optimization=false, incremental_redundancy=false)
    result2 = AS(rpd; tol=1e-9, use_optimization=true, incremental_redundancy=true)
    
    println(@sprintf("Original vertices:  %d", size(result1, 1)))
    println(@sprintf("Optimized vertices: %d", size(result2, 1)))
    
    # Expected result
    expected = [3.0 3.0; 3.0 9.75; 9.0 9.0; 9.75 3.0]
    
    function check_result(vertices, name)
        matches = 0
        for i in 1:size(expected, 1)
            for j in 1:size(vertices, 1)
                if abs(vertices[j,1] - expected[i,1]) < 1e-6 && abs(vertices[j,2] - expected[i,2]) < 1e-6
                    matches += 1
                    break
                end
            end
        end
        correct = matches == size(expected, 1)
        println(@sprintf("%s: %s (%d/%d expected points)", name, correct ? "✅" : "❌", matches, size(expected, 1)))
        return correct
    end
    
    correct1 = check_result(result1, "Original")
    correct2 = check_result(result2, "Optimized")
    
    if correct1 && correct2
        println("✅ Both methods produce correct results!")
    else
        println("❌ Correctness issue detected!")
    end
    
    return (original_time, optimized_time, speedup)
end

# Run test
try
    run_benchmark()
catch e
    println("Error: $e")
    # Show stack trace for debugging
    Base.show_backtrace(stdout, catch_backtrace())
end