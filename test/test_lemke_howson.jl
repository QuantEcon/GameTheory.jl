@testset "lemke_howson.jl" begin

    function test_lh_case(case)
        g = case[:g]
        NEs_dict = case[:NEs_dict]
        converged = case[:converged]

        for k in keys(NEs_dict)
            NE_compupted, res =
                lemke_howson(g, init_pivot=k, full_output=Val(true))
            @test isapprox_act_profs(NE_compupted, NEs_dict[k])
            @test res.converged == converged
        end
    end

    @testset "Basic" begin
        # From von Stengel 2007 in Algorithmic Game Theory
        A = [3 3; 2 5; 0 6]
        B = [3 2 3; 2 6 1]
        g = NormalFormGame(Player(A), Player(B))
        NEs_dict = Dict(
            1 => ([1., 0., 0.], [1., 0.]),  # init_pivot => NE
            2 => ([0., 1/3, 2/3], [1/3, 2/3])
        )

        NE = @inferred lemke_howson(g)
        @test is_nash(g, NE)
        NE = @inferred lemke_howson(g, init_pivot=5)
        @test is_nash(g, NE)

        NE, res = @inferred lemke_howson(g, full_output=Val(true))
        @test is_nash(g, res.NE)
        @test res.converged

        T1 = Rational{BigInt}
        T2 = BigFloat
        g_1 = NormalFormGame(T1, g)
        NE_1 = @inferred lemke_howson(g_1)
        @test is_nash(g_1, NE_1)
        @test eltype(NE_1[1]) == T2

        case = (g = g, NEs_dict = NEs_dict, converged = true)
        test_lh_case(case)
    end

    @testset "Degenerate games" begin
        cases = []

        # From von Stengel 2007 in Algorithmic Game Theory
        A = [3 3; 2 5; 0 6]
        B = [3 2 3; 3 6 1]
        push!(cases, (
            g = NormalFormGame(Player(A), Player(B)),
            NEs_dict = Dict(1 => ([0., 1/3, 2/3], [1/3, 2/3])),
            converged = true
        ))

        # == Examples of cycles by "ad hoc" tie breaking rules == #

        # Example where tie breaking that picks the variable with
        # the smallest row index in the tableau leads to cycling
        A = [0 0 0;
             0 1 1;
             1 1 0]
        B = [1 0 1;
             1 1 0;
             0 0 2]
        push!(cases, (
            g = NormalFormGame(Player(A), Player(B)),
            NEs_dict = Dict(1 => ([0., 2/3, 1/3], [0., 1., 0.])),
            converged = true
        ))

        # Example where tie breaking that picks the variable with
        # the smallest variable index in the tableau leads to cycling
        perm = [3, 1, 2]
        C = A[:, perm]
        D = B[perm, :]
        push!(cases, (
            g = NormalFormGame(Player(C), Player(D)),
            NEs_dict = Dict(1 => ([0., 2/3, 1/3], [0., 0., 1.])),
            converged = true
        ))

        test_lh_case.(cases)
    end

    @testset "Capping" begin
        A = [3 3; 2 5; 0 6]
        B = [3 2 3; 2 6 1]
        g = NormalFormGame(Player(A), Player(B))
        m, n = g.nums_actions
        max_iter = 10^6  # big number

        for k in 1:m+n
            (NE1, res1) = @inferred(
                lemke_howson(g; init_pivot=k, max_iter=max_iter,
                             capping=nothing, full_output = Val(true))
            )
            (NE2, res2) = @inferred(
                lemke_howson(g; init_pivot=k, max_iter = max_iter,
                             capping=max_iter, full_output = Val(true))
            )
            @test isapprox(NE1[1], NE2[1])
            @test isapprox(NE1[2], NE2[2])
            @test res1.init == res2.init
        end

        init_pivot = 2
        max_iter = m+n
        (NE, res) = lemke_howson(g; init_pivot=init_pivot, max_iter=max_iter,
                                 capping = 1, full_output=Val(true))
        @test res.num_iter == max_iter
        @test res.init == init_pivot-1
    end

    @testset "Invalid init_pivot" begin
        A = [3 3; 2 5; 0 6]
        B = [3 2 3; 2 6 1]
        g = NormalFormGame(Player(A), Player(B))
        m, n = g.nums_actions

        @test_throws ArgumentError lemke_howson(g; init_pivot=-1)
        @test_throws ArgumentError lemke_howson(g; init_pivot=0)
        @test_throws ArgumentError lemke_howson(g; init_pivot=m+n+1)
    end

end
