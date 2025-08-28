#!/usr/bin/env julia

using Pkg
Pkg.activate(".")
using GameTheory
using Printf

"""
Comprehensive Performance and Correctness Analysis 
Comparing AS algorithm versions 694f4e5 vs dbbb640
"""

function comprehensive_benchmark()
    println("🔬 AS Algorithm Performance Analysis")
    println("=" ^ 70)
    println("Comparing versions:")
    println("  • 694f4e5 (current): With optimizations from PR #65 comment 5")
    println("  • dbbb640 (previous): Original implementation")
    println("  • Using optimization flags to simulate version differences")
    println()
    
    # Test case
    pd_payoff = [9.0 1.0; 10.0 3.0]
    nfg = NormalFormGame(Player(pd_payoff), Player(pd_payoff))
    rpd = RepeatedGame(nfg, 0.75)
    
    println("Test Case: 2x2 Prisoner's Dilemma with δ=0.75")
    println("Expected vertices: [3.0, 3.0], [3.0, 9.75], [9.0, 9.0], [9.75, 3.0]")
    println()
    
    # Performance measurement
    function time_method(name, use_opt, incr_red, n_runs=5)
        times = Float64[]
        result = nothing
        
        print("$name: ")
        for i in 1:n_runs
            GC.gc()
            t = @elapsed result = AS(rpd; tol=1e-9, 
                                   use_optimization=use_opt, 
                                   incremental_redundancy=incr_red)
            push!(times, t)
            print(".")
        end
        
        avg_time = sum(times) / length(times)
        min_time = minimum(times)
        
        println(@sprintf(" %.4fs (avg), %.4fs (min)", avg_time, min_time))
        return min_time, result
    end
    
    # Test both methods
    println("📊 PERFORMANCE MEASUREMENT (5 runs each):")
    println("-" ^ 50)
    original_time, original_result = time_method(
        "Original (dbbb640 equiv.)", false, false)
    
    optimized_time, optimized_result = time_method(
        "Optimized (694f4e5)      ", true, true)
    
    speedup = original_time / optimized_time
    
    println()
    println("⚡ PERFORMANCE SUMMARY:")
    println("-" ^ 30)
    println(@sprintf("Original time:   %.4f seconds", original_time))
    println(@sprintf("Optimized time:  %.4f seconds", optimized_time))
    println(@sprintf("Speedup:         %.2fx", speedup))
    
    if speedup > 2.0
        println("🚀 Excellent speedup - optimizations very effective!")
    elif speedup > 1.2
        println("✅ Good speedup - optimizations effective")
    elseif speedup > 1.0
        println("✅ Modest speedup - optimizations helpful")
    else
        println("📉 Optimizations slower than original")
    end
    
    # Correctness check
    println()
    println("🔍 CORRECTNESS ANALYSIS:")
    println("-" ^ 30)
    
    expected = [3.0 3.0; 3.0 9.75; 9.0 9.0; 9.75 3.0]
    
    function analyze_correctness(result, name)
        println("$name results:")
        println("  Size: $(size(result))")
        println("  Vertices:")
        for i in 1:size(result, 1)
            println(@sprintf("    [%.4f, %.4f]", result[i,1], result[i,2]))
        end
        
        # Check matches
        matches = 0
        for i in 1:size(expected, 1)
            found = false
            for j in 1:size(result, 1)
                if abs(result[j,1] - expected[i,1]) < 1e-6 && abs(result[j,2] - expected[i,2]) < 1e-6
                    found = true
                    matches += 1
                    break
                end
            end
        end
        
        correct = matches == size(expected, 1)
        println(@sprintf("  Correctness: %s (%d/%d expected vertices found)", 
                correct ? "✅ CORRECT" : "❌ INCORRECT", matches, size(expected, 1)))
        println()
        return correct
    end
    
    original_correct = analyze_correctness(original_result, "Original")
    optimized_correct = analyze_correctness(optimized_result, "Optimized")
    
    # Overall assessment
    println("🎯 OVERALL ASSESSMENT:")
    println("=" ^ 40)
    
    if original_correct && optimized_correct
        println("✅ Correctness: Both methods produce correct results")
        if speedup > 1.1
            println("🎉 CONCLUSION: Optimizations successful!")
            println("   The optimizations from PR #65 provide significant benefit")
            println("   without sacrificing correctness.")
        else
            println("📊 CONCLUSION: Optimizations provide minimal benefit")
            println("   But they don't hurt correctness.")
        end
    else
        println("❌ CORRECTNESS ISSUE DETECTED!")
        
        if !original_correct
            println("   Original method has correctness issues!")
        end
        
        if !optimized_correct
            println("   ⚠️  CRITICAL: Optimized method produces wrong results!")
            println("   The optimizations introduce bugs that need to be fixed.")
            println("   Performance improvement is meaningless if results are wrong.")
        end
        
        println("🔧 RECOMMENDATION: Fix correctness issues before considering performance")
    end
    
    # Test individual optimizations if there's a correctness issue
    if !optimized_correct && original_correct
        println()
        println("🔧 DEBUGGING: Testing individual optimizations...")
        println("-" ^ 50)
        
        _, opt1_result = time_method("use_optimization only     ", true, false, 3)
        _, opt2_result = time_method("incremental_redundancy only", false, true, 3)
        
        opt1_correct = analyze_correctness(opt1_result, "use_optimization only")
        opt2_correct = analyze_correctness(opt2_result, "incremental_redundancy only")
        
        if opt1_correct && !opt2_correct
            println("🎯 Root cause: incremental_redundancy=true has bugs")
        elif !opt1_correct && opt2_correct
            println("🎯 Root cause: use_optimization=true has bugs")  
        elif !opt1_correct && !opt2_correct
            println("🎯 Root cause: Both optimizations have issues")
        else
            println("🤔 Root cause: Issues only occur when both optimizations are combined")
        end
    end
    
    return (
        original_time=original_time,
        optimized_time=optimized_time,
        speedup=speedup,
        original_correct=original_correct,
        optimized_correct=optimized_correct
    )
end

# Run the analysis
try
    results = comprehensive_benchmark()
    println("\n" * "=" * 70)
    println("Analysis completed successfully!")
catch e
    println("Error during analysis: $e")
    Base.show_backtrace(stdout, catch_backtrace())
end