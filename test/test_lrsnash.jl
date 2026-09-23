@testset "lrsnash.jl" begin

    @testset "test 3 by 2 non-degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([3 3; 2 5; 0 6]),
                           Player([3 2 3; 2 6 1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game(Rational)" begin
        g = NormalFormGame(Player([3//1 3//1; 2//1 5//1; 0//1 6//1]),
                           Player([3//1 2//1 3//1; 2//1 6//1 1//1]))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test 3 by 2 non-degenerate normal form game($T)" for T in
            [BigInt, Int32, UInt8, Rational{UInt8}]
        g = NormalFormGame(Player(T.([3 3; 2 5; 0 6])),
                           Player(T.([3 2 3; 2 6 1])))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test payoffs beyond Int64" begin
        c = big(2)^70
        g = NormalFormGame(Player(c .* [3 3; 2 5; 0 6] .- c),
                           Player(c .* [3 2 3; 2 6 1] .+ 1))
        NEs = [([1//1, 0//1, 0//1], [1//1, 0//1]),
               ([4//5, 1//5, 0//1], [2//3, 1//3]),
               ([0//1, 1//3, 2//3], [1//3, 2//3])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test 2 by 2 coordination game(Bool)" begin
        g = NormalFormGame(Player([true false; false true]),
                           Player([true false; false true]))
        NEs = [([1//1, 0//1], [1//1, 0//1]),
               ([1//2, 1//2], [1//2, 1//2]),
               ([0//1, 1//1], [0//1, 1//1])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "test Float64 payoffs rejected" begin
        g = NormalFormGame(Player([3. 3.; 2. 5.; 0. 6.]),
                           Player([3. 2. 3.; 2. 6. 1.]))
        @test_throws MethodError lrsnash(g)
    end

    @testset "test 3 by 2 degenerate normal form game(Int)" begin
        g = NormalFormGame(Player([1 -1; -1 1; 0 0]),
                           Player([1 0 0; 0 0 0]))
        NEs = [([1, 0, 0], [1, 0]),
               ([0, 1, 0], [0, 1]),
               ([0, 1, 0], [1//2, 1//2]),
               ([0, 0, 1], [1//2, 1//2])]
        NEs_computed = @inferred(lrsnash(g))

        @test sort(NEs_computed) == sort(NEs)
    end

end
