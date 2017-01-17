@testset "Testing Repeated Game functionality" begin

    pd_payoff = [9.0 1.0
                 10.0 3.0]

    A = Player(pd_payoff)
    B = Player(pd_payoff)
    nfg = NormalFormGame((A, B))

    # Tests construction of repeated game
    rpd = RepeatedGame(nfg, 0.75)
    C, H, Z = Games.initialize_sg_hpl(4, [0.0, 0.0], 1.0)

    #
    # Test various helper functions
    #
    @testset "Testing flow utility computations" begin
        @test abs(flow_u_1(rpd, 1, 1) - 9.0) < 1e-14
        @test abs(flow_u_2(rpd, 2, 2) - 3.0) < 1e-14
        @test maxabs(flow_u(rpd, 2, 2) - [3.0, 3.0]) < 1e-14
    end

    @testset "Testing best deviation computations" begin
        @test best_dev_i(rpd, 1, 1) == 2
    end

    @testset "Testing unit circle function" begin
        H = Games.unitcircle(4)
        points = [1.0 0.0
                  0.0 1.0
                  -1.0 0.0
                  0.0 -1.0]

        @test maxabs(H - points) < 1e-12
    end

    @testset "Testing subgradient and hyperplane level initialize" begin
        C, H, Z = Games.initialize_sg_hpl(4, [0.0, 0.0], 1.0)

        @test maxabs(C - ones(4)) < 1e-12
        @test maxabs(H - Z') < 1e-12
    end

    @testset "Testing worst value computation" begin
        @test abs(worst_value_i(rpd, H, C, 1) + 1.0) < 1e-12
    end

    #
    # Test the actual computation
    #
    @testset "Testing outer approximation" begin
        vertices = outerapproximation(rpd; nH=128, maxiter=150, tol=1e-9)
        p_in_v = [vertices[i, :] for i in 1:size(vertices, 1)]

        mybools = [all(isapprox.([3.0, 3.0], p)) for p in p_in_v]
        @test any(mybools)
    end

end

