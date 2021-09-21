using CDDLib

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
        dict = Dict(:nH=>64, :maxiter=>150, :tol=>1e-9)
        plib = CDDLib.Library()
        for kwargs in [dict, merge(dict, Dict(:plib=>plib))]
            vertices = @inferred(outerapproximation(rpd; kwargs...))
            p_in_v = [vertices[i, :] for i in 1:size(vertices, 1)]

            mybools = [all(isapprox.([3.0, 3.0], p)) for p in p_in_v]
            @test any(mybools)
        end
    end

end

