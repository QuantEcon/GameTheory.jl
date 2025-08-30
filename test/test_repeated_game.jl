@testset "Testing Repeated Game functionality" begin

    pd_payoff = [9.0 1.0
                 10.0 3.0]

    A = Player(pd_payoff)
    B = Player(pd_payoff)
    nfg = NormalFormGame((A, B))

    # Tests construction of repeated game
    rpd = RepeatedGame(nfg, 0.75)
    C, H, Z = GameTheory.initialize_sg_hpl(4, [0.0, 0.0], 1.0)

    #
    # Test various helper functions
    #
    @testset "Testing flow utility computations" begin
        @test abs(flow_u_1(rpd, 1, 1) - 9.0) < 1e-14
        @test abs(flow_u_2(rpd, 2, 2) - 3.0) < 1e-14
        @test maximum(abs, flow_u(rpd, 2, 2) - [3.0, 3.0]) < 1e-14
    end

    @testset "Testing best deviation computations" begin
        @test best_dev_i(rpd, 1, 1) == 2
    end

    @testset "Testing unit circle function" begin
        H = GameTheory.unitcircle(4)
        points = [1.0 0.0
                  0.0 1.0
                  -1.0 0.0
                  0.0 -1.0]

        @test maximum(abs, H - points) < 1e-12
    end

    @testset "Testing subgradient and hyperplane level initialize" begin
        C, H, Z = GameTheory.initialize_sg_hpl(4, [0.0, 0.0], 1.0)

        @test maximum(abs, C - ones(4)) < 1e-12
        @test maximum(abs, H - Z') < 1e-12
    end

    @testset "Testing worst value computation" begin
        @test abs(worst_value_i(rpd, H, C, 1) + 1.0) < 1e-12
    end

    #
    # Test the actual computation
    #
    @testset "Testing outer approximation" begin
        kwargs = Dict(:nH=>64, :maxiter=>150, :tol=>1e-9)
        vertices = @inferred(outerapproximation(rpd; kwargs...))
        p_in_v = [vertices[i, :] for i in 1:size(vertices, 1)]

        mybools = [all(isapprox.([3.0, 3.0], p)) for p in p_in_v]
        @test any(mybools)
    end

    #
    # Test AS algorithm
    #
    @testset "Testing AS algorithm" begin
        # Helper function to test if each row of vertices approximately matches some row in expected
        function vertices_match_expected(vertices, expected; tol=1e-4)
            # Convert to Float64 for comparison if needed
            vertices_float = eltype(vertices) <: AbstractFloat ? vertices : Float64.(vertices)
            expected_float = eltype(expected) <: AbstractFloat ? expected : Float64.(expected)
            
            # Check that each row in vertices is approximately equal to some row in expected
            for i in 1:size(vertices_float, 1)
                found_match = false
                for j in 1:size(expected_float, 1)
                    if maximum(abs, vertices_float[i, :] - expected_float[j, :]) < tol
                        found_match = true
                        break
                    end
                end
                if !found_match
                    return false
                end
            end
            return true
        end
        
        vertices = @inferred(AS(rpd; tol=1e-9))

        pts_sorted = [3.0 3.0;
                      3.0 9.75;
                      9.0 9.0;
                      9.75 3.0]
        
        # Test that each row of vertices is approximately equal to some row of pts_sorted
        @test vertices_match_expected(vertices, pts_sorted)

        @testset "AS with Int payoffs" begin
            nfg_int = NormalFormGame(Int, nfg)
            rpd_int = RepeatedGame(nfg_int, 0.75)
            vertices = @inferred(AS(rpd_int; tol=1e-9))
            @test vertices_match_expected(vertices, pts_sorted)

            vertices_u = @inferred(AS(rpd_int; tol=1e-9, u=[0, 0]))
            @test vertices_match_expected(vertices_u, pts_sorted)
        end

        @testset "AS with Rational payoffs" begin
            nfg_rat = NormalFormGame(Rational{Int}, nfg)
            rpd_rat = RepeatedGame(nfg_rat, 0.75)
            vertices = @inferred(AS(rpd_rat; tol=1e-9))
            @test vertices_match_expected(vertices, pts_sorted)
        end

        @testset "AS with rational payoffs and rational delta" begin
            # Test the case described in the issue with both rational payoffs and delta
            nfg_rat = NormalFormGame(Rational{Int}, nfg)
            rpd_rat = RepeatedGame(nfg_rat, 3//4)  # Rational delta
            vertices = @inferred(AS(rpd_rat; tol=1e-9))
            @test vertices_match_expected(vertices, pts_sorted)
            # For rational case, the result should be Matrix{Rational{BigInt}}
            @test eltype(vertices) == Rational{BigInt}
        end

        @testset "AS with Int payoffs and rational delta" begin
            # Test the case with Int payoffs but rational delta (should use exact arithmetic)
            nfg_int = NormalFormGame(Int, nfg)
            rpd_int_rat = RepeatedGame(nfg_int, 3//4)  # Rational delta with Int payoffs
            vertices = @inferred(AS(rpd_int_rat; tol=1e-9))
            @test vertices_match_expected(vertices, pts_sorted)
            # Should also use exact arithmetic and return Rational{BigInt}
            @test eltype(vertices) == Rational{BigInt}
        end

        @testset "AS with verbose output" begin
            # Test verbose parameter
            rpd_test = RepeatedGame(nfg, 0.75)
            # Should run without error and print convergence message
            vertices = @inferred(AS(rpd_test; tol=1e-9, verbose=true))
            @test size(vertices) == size(pts_sorted)
        end

        @testset "uniquetolrows function" begin
            # Test the uniquetolrows utility function 
            V = [1.0001 2.0002; 1.0 2.0; 3.0 4.0; 1.00009 2.00008]
            tol = 1e-3
            V_unique = uniquetolrows(V, tol)
            # Should remove the near-duplicate rows (rows 1, 2, 4 are all within tolerance)
            @test size(V_unique, 1) == 2  # Two duplicates should be removed
            
            # Test with the example from the issue
            tol = 1e-9
            rpd_test = RepeatedGame(nfg, 0.75)
            V_test = AS(rpd_test; tol=tol)
            V_approx = uniquetolrows(V_test, tol)
            @test size(V_approx, 1) <= size(V_test, 1)  # Should not increase size
        end
    end

end
