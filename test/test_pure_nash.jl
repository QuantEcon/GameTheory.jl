using Combinatorics


@testset "Testing Pure Nash Equilibrium Routines" begin

    # Pure Action Nash equilibrium
    @testset "2x2 Game with 1 Pure Action Nash equilibrium" begin
        A = [9.0 1.0
             10.0 3.0]

        nfg = NormalFormGame(A)
        ne = pure_nash(nfg)

        @test ne == [(2, 2)]
    end

    @testset "Matching Pennies Game with 0 Pure Action Nash equilibrium" begin
        MP = [1.0 -1.0
              -1.0 1.0]
        p1 = Player(MP)
        p2 = Player(-MP)

        g_MP = NormalFormGame(p1, p2)
        ne = pure_nash(g_MP)

        @test ne == []
    end

    @testset "Coordination game with 2 Pure Action Nash equilibria" begin
        Coo = [4.0 0.0
               3.0 2.0]

        g_Coo = NormalFormGame(Coo)
        ne = pure_nash(g_Coo)

        @test sort(ne) == sort([(1,1); (2,2)])
    end

    @testset "Coordination game with 2 Pure Action Nash equilibria but only find 1" begin
        Coo = [4.0 0.0
               3.0 2.0]

        g_Coo = NormalFormGame(Coo)
        ne = pure_nash(g_Coo; ntofind=1)

        @test length(ne) == 1
    end

    @testset "Unanimity Game with more than two players" begin
        N = 4
        a, b = 1, 2
        g_Unanimity = NormalFormGame(tuple(fill(2, N)...))
        g_Unanimity[fill(1, N)...] = fill(a, N)
        g_Unanimity[fill(2, N)...] = fill(b, N)

        Unanimity_NE = [tuple(fill(1, N)...)]
        for k in 2:N-2
            for ind in combinations(1:N, k)
                a = fill(1, N)
                a[ind] = 2
                push!(Unanimity_NE, tuple(a...))
            end
        end
        push!(Unanimity_NE, tuple(fill(2, N)...))

        ne = pure_nash(g_Unanimity)

        @test sort(ne) == sort(Unanimity_NE)
    end

    @testset "Tolerance" begin
        epsilon = 1e-08

        g = NormalFormGame((2, 2))
        g[1, 1] = [1, 1]
        g[1, 2] = [-2, 1 + epsilon]
        g[2, 1] = [1 + epsilon, -2]
        g[2, 2] = [0, 0];

        NEs = [[(2, 2)]]
        epsilon_NEs = [[(1, 1); (2, 2)]]

        for (tol, answer) in zip([0 epsilon], [NEs epsilon_NEs])
            @test sort(pure_nash(g, tol=tol)) ==  sort(answer)
        end
    end

    @testset "Trivial game with 1 player" begin
        n = 3
        g1 = NormalFormGame(Player(collect(1:n)))
        @test pure_nash(g1) == [(n,)]
        @test sort(pure_nash(g1, tol=1.)) == [(n-1,), (n,)]
    end

end
