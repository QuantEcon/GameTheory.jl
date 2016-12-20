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

end
