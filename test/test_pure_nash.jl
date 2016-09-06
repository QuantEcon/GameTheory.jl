@testset "Testing Pure Nash Equilibrium Routines" begin

    # Pure strategy Nash equilibrium
    @testset "2x2 Pure strategy Nash equilibrium" begin
        A = [9.0 1.0
             10.0 3.0]
        p1 = Player(A)
        p2 = Player(A)

        nfg = NormalFormGame(p1, p2)
        psne = pure_nash(nfg)

        @test psne == [(2, 2)]
    end

end

