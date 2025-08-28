#!/usr/bin/env julia

"""
Simple performance comparison for AS algorithm optimizations from PR #65.

Compares timing and results between:
- Original method (use_optimization=false, incremental_redundancy=false) 
- Both optimizations enabled (use_optimization=true, incremental_redundancy=true)
"""

using LinearAlgebra
using Printf

# We need to be in the right directory to import GameTheory
cd("/home/runner/work/GameTheory.jl/GameTheory.jl")

# Manual import since we're having package issues
include("src/GameTheory.jl")
using .GameTheory

function create_test_games()
    games = []
    
    # 2x2 Prisoner's Dilemma (from tests)
    pd_payoff = [9.0 1.0; 10.0 3.0]
    A = Player(pd_payoff)
    B = Player(pd_payoff)
    nfg_pd = NormalFormGame((A, B))
    push!(games, ("2x2 PD Float64", RepeatedGame(nfg_pd, 0.75)))
    push!(games, ("2x2 PD Rational", RepeatedGame(NormalFormGame(Rational{Int}, nfg_pd), 3//4)))
    
    # 3x3 game for more complex testing
    payoff_3x3 = [5.0 2.0 1.0; 6.0 4.0 2.0; 3.0 5.0 4.0]
    A3 = Player(payoff_3x3)
    B3 = Player(payoff_3x3')
    nfg_3x3 = NormalFormGame((A3, B3))
    push!(games, ("3x3 Game Float64", RepeatedGame(nfg_3x3, 0.8)))
    
    return games
end

function time_as_variant(rpd, name, use_opt, incr_red; maxiter=100, n_runs=3)
    """Run AS with specific optimization settings and measure performance"""
    
    println("  Testing $name...")
    
    times = Float64[]
    vertices = 0
    
    for i = 1:n_runs
        try
            GC.gc()  # Clean up memory
            start_time = time_ns()
            result = AS(rpd; maxiter=maxiter, tol=1e-7, 
                       use_optimization=use_opt, 
                       incremental_redundancy=incr_red,
                       verbose=false)
            end_time = time_ns()
            
            elapsed = (end_time - start_time) / 1e9
            push!(times, elapsed)
            vertices = size(result, 1)
            
        catch e
            println("    ❌ Failed on run $i: $e")
            return (time=Inf, vertices=0, success=false)
        end
    end
    
    avg_time = sum(times) / length(times)
    min_time = minimum(times)
    
    println(@sprintf "    ⏱  Avg Time: %.4f s (min: %.4f s)", avg_time, min_time)
    println(@sprintf "    📊 Vertices: %d", vertices)
    
    return (time=avg_time, min_time=min_time, vertices=vertices, success=true)
end

function run_performance_comparison()
    println("🔬 AS Algorithm Performance Comparison")
    println("=" ^ 60)
    println("Comparing commit 694f4e5 (with optimizations) vs dbbb640 (original)")
    println("Using optimization flags to simulate the difference")
    println()
    
    games = create_test_games()
    
    # Test configurations - simulate old vs new behavior
    configs = [
        ("Original Method", false, false),       # Simulates dbbb640 behavior  
        ("With Optimizations", true, true)      # Current 694f4e5 behavior
    ]
    
    results = Dict()
    
    for (game_name, rpd) in games
        println("\n🎮 Game: $game_name")
        println("-" ^ 40)
        
        game_results = Dict()
        
        for (config_name, use_opt, incr_red) in configs
            result = time_as_variant(rpd, config_name, use_opt, incr_red)
            game_results[config_name] = result
        end
        
        results[game_name] = game_results
        
        # Print comparison for this game
        original_result = game_results["Original Method"]
        optimized_result = game_results["With Optimizations"]
        
        if original_result.success && optimized_result.success
            speedup = original_result.time / optimized_result.time
            println()
            if speedup > 1.0
                println(@sprintf "  📈 Performance Improvement: %.2fx faster with optimizations", speedup)
            elseif speedup < 1.0
                println(@sprintf "  📉 Performance Regression: %.2fx slower with optimizations", 1/speedup)
            else
                println("  ➡️  Similar performance")
            end
            
            println(@sprintf "  📊 Vertex comparison: %d (original) vs %d (optimized)", 
                    original_result.vertices, optimized_result.vertices)
        end
        println()
    end
    
    # Overall summary
    println("\n" * "=" * 60)
    println("📋 OVERALL PERFORMANCE SUMMARY")  
    println("=" * 60)
    
    all_speedups = Float64[]
    
    for (game_name, game_results) in results
        original_result = game_results["Original Method"]
        optimized_result = game_results["With Optimizations"]
        
        if original_result.success && optimized_result.success
            speedup = original_result.time / optimized_result.time
            push!(all_speedups, speedup)
            
            status = speedup > 1.05 ? "✅ FASTER" : speedup < 0.95 ? "📉 SLOWER" : "➡️ SIMILAR"
            println(@sprintf "%s: %s (%.2fx)", game_name, status, 
                    speedup > 1.0 ? speedup : 1/speedup)
        else
            println("$game_name: ⚠ Incomplete benchmark")
        end
    end
    
    if !isempty(all_speedups)
        avg_speedup = sum(all_speedups) / length(all_speedups)
        println()
        println(@sprintf "Average Performance Change: %.2fx", avg_speedup)
        
        if avg_speedup > 1.05
            println("🎉 Overall: Optimizations provide meaningful speed improvement!")
        elseif avg_speedup < 0.95
            println("⚠️  Overall: Optimizations appear to slow down performance")
        else
            println("📊 Overall: Optimizations have minimal performance impact")
        end
    end
    
    return results
end

# Run the benchmark
if abspath(PROGRAM_FILE) == @__FILE__
    run_performance_comparison()
end