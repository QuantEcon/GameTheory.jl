#!/usr/bin/env julia

"""
Performance Benchmark for AS Algorithm Optimizations

This script compares the performance between:
- Version dbbb640: Original implementation (use_optimization=false, incremental_redundancy=false)  
- Version 694f4e5: Optimized implementation (use_optimization=true, incremental_redundancy=true)

The optimizations implement suggestions from PR #65 comment 5:
1. Optimization-based IC boundary intersection (use_optimization=true)
2. Incremental redundancy removal (incremental_redundancy=true)
"""

using GameTheory
using BenchmarkTools
using Printf

# Test configurations
struct TestCase
    name::String
    nfg::NormalFormGame
    delta::Union{Float64, Rational{Int}}
end

function create_test_cases()
    cases = TestCase[]
    
    # Small 2x2 Prisoner's Dilemma (Float64)
    pd_payoff = [9.0 1.0; 10.0 3.0]
    nfg_2x2_float = NormalFormGame(Player(pd_payoff), Player(pd_payoff))
    push!(cases, TestCase("2x2 Float64 δ=0.75", nfg_2x2_float, 0.75))
    
    # Small 2x2 with Int payoffs (Float64 delta)
    nfg_2x2_int = NormalFormGame(Int, nfg_2x2_float)
    push!(cases, TestCase("2x2 Int δ=0.75", nfg_2x2_int, 0.75))
    
    # Small 2x2 with Rational payoffs and rational delta
    nfg_2x2_rat = NormalFormGame(Rational{Int}, nfg_2x2_float)
    push!(cases, TestCase("2x2 Rational δ=3//4", nfg_2x2_rat, 3//4))
    
    # Larger 3x3 game (Float64)
    payoff_3x3 = [8.0 2.0 1.0; 9.0 6.0 3.0; 4.0 7.0 5.0]
    nfg_3x3_float = NormalFormGame(Player(payoff_3x3), Player(payoff_3x3))
    push!(cases, TestCase("3x3 Float64 δ=0.8", nfg_3x3_float, 0.8))
    
    # Larger 3x3 with rational delta (more demanding)
    push!(cases, TestCase("3x3 Float64 δ=4//5", nfg_3x3_float, 4//5))
    
    return cases
end

function benchmark_as_version(test_case::TestCase, use_optimization::Bool, incremental_redundancy::Bool; samples=3)
    rpd = RepeatedGame(test_case.nfg, test_case.delta)
    
    # Run benchmark
    result = @benchmark AS($rpd; 
                           tol=1e-9, 
                           use_optimization=$use_optimization,
                           incremental_redundancy=$incremental_redundancy) samples=samples evals=1
    
    return result
end

function run_performance_comparison()
    println("=" ^ 80)
    println("AS Algorithm Performance Comparison")
    println("=" ^ 80)
    println("Comparing optimizations from PR #65 comment 5:")
    println("  - Original (dbbb640): use_optimization=false, incremental_redundancy=false")
    println("  - Optimized (694f4e5): use_optimization=true, incremental_redundancy=true")
    println()
    
    test_cases = create_test_cases()
    
    # Results storage
    results = []
    
    for (i, test_case) in enumerate(test_cases)
        println("[$i/$(length(test_cases))] Testing: $(test_case.name)")
        println("-" ^ 60)
        
        try
            # Benchmark original version (dbbb640 equivalent)
            println("  Running original version...")
            original_result = benchmark_as_version(test_case, false, false; samples=3)
            original_time = minimum(original_result.times) / 1e9  # Convert to seconds
            
            # Benchmark optimized version (694f4e5)
            println("  Running optimized version...")
            optimized_result = benchmark_as_version(test_case, true, true; samples=3)
            optimized_time = minimum(optimized_result.times) / 1e9  # Convert to seconds
            
            # Calculate speedup
            speedup = original_time / optimized_time
            
            # Store results
            push!(results, (
                name = test_case.name,
                original_time = original_time,
                optimized_time = optimized_time,
                speedup = speedup
            ))
            
            # Print results
            println("  Original time:  $(@sprintf("%.4f", original_time)) seconds")
            println("  Optimized time: $(@sprintf("%.4f", optimized_time)) seconds")
            println("  Speedup:        $(@sprintf("%.2f", speedup))x")
            
            if speedup > 1.0
                println("  ✅ Optimization successful - $(@sprintf("%.1f", (speedup-1)*100))% faster")
            elseif speedup < 0.9
                println("  ⚠️  Optimization slower - $(@sprintf("%.1f", (1-speedup)*100))% slower")
            else
                println("  ≈  Similar performance")
            end
            
        catch e
            println("  ❌ Error: $e")
            push!(results, (
                name = test_case.name,
                original_time = NaN,
                optimized_time = NaN,
                speedup = NaN
            ))
        end
        
        println()
    end
    
    # Summary
    println("=" ^ 80)
    println("PERFORMANCE SUMMARY")
    println("=" ^ 80)
    
    valid_results = filter(r -> !isnan(r.speedup), results)
    
    if !isempty(valid_results)
        avg_speedup = sum(r.speedup for r in valid_results) / length(valid_results)
        max_speedup = maximum(r.speedup for r in valid_results)
        min_speedup = minimum(r.speedup for r in valid_results)
        
        println(@sprintf("Average speedup: %.2fx", avg_speedup))
        println(@sprintf("Best speedup:    %.2fx (%s)", max_speedup, 
                valid_results[argmax([r.speedup for r in valid_results])].name))
        println(@sprintf("Worst speedup:   %.2fx (%s)", min_speedup, 
                valid_results[argmin([r.speedup for r in valid_results])].name))
        
        println()
        println("Detailed Results:")
        println(@sprintf("%-25s %12s %12s %8s", "Test Case", "Original", "Optimized", "Speedup"))
        println("-" ^ 60)
        for result in results
            if !isnan(result.speedup)
                println(@sprintf("%-25s %8.4fs %8.4fs %7.2fx", 
                        result.name, result.original_time, result.optimized_time, result.speedup))
            else
                println(@sprintf("%-25s %12s %12s %8s", result.name, "ERROR", "ERROR", "N/A"))
            end
        end
        
        println()
        if avg_speedup > 1.1
            println("🚀 Overall conclusion: Optimizations provide significant performance improvement!")
            println("   The optimizations from PR #65 comment 5 are effective.")
        elseif avg_speedup > 1.0
            println("✅ Overall conclusion: Optimizations provide modest performance improvement.")
        else
            println("⚠️  Overall conclusion: Optimizations do not improve performance on average.")
            println("   May need further investigation or different test cases.")
        end
    else
        println("❌ No valid benchmarks completed. All tests failed.")
    end
    
    println()
    println("Note: Times shown are minimum times from 3 benchmark samples.")
    println("      Speedup = original_time / optimized_time")
    println("=" ^ 80)
    
    return results
end

# Run the benchmark if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_performance_comparison()
end