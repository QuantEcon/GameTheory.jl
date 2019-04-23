# ------------------------- #
# Testing local interaction #
# ------------------------- #


@testset "Testing localint.jl" begin
    payoff_matrix = [4 0; 2 3]
    adj_matrix = [0 1 3; 2 0 1; 3 2 0]
    init_actions = (1, 1, 2)

    @testset "LocalInteraction from NormalFormGame" begin
        game = NormalFormGame(payoff_matrix)
        li = LocalInteraction(game, adj_matrix)

        @test li.players == ntuple(i -> Player(payoff_matrix), 3)
        @test li.num_actions == 2
    end

    @testset "Testing local interaction with simultaneous revision" begin
        li = LocalInteraction(payoff_matrix, adj_matrix)

        @test @inferred(play(li, init_actions) == (2, 1, 1))
        @test @inferred(play(li, init_actions, num_reps=2) == (1, 2, 2))
        @test @inferred(time_series(li, 3, init_actions) == [1 2 1; 1 1 2; 2 1 2])
        @test all(time_series(li, 3) .<= li.num_actions)
    end

    @testset "Testing local interaction with asynchronous revision" begin
        seed = 1234
        li = LocalInteraction(payoff_matrix, adj_matrix, AsynchronousRevision())

        @test @inferred(play(li, init_actions, 1) == (2, 1, 2))
        @test @inferred(play(li, init_actions, [1,2]) == (2, 1, 2))
        @test @inferred(play(li, init_actions, [1,2], num_reps=2) == (2, 2, 2))
        @test all(time_series(MersenneTwister(seed), li, 3, init_actions)
                  .<= li.num_actions)
        @test all(time_series(MersenneTwister(seed), li, 3, [1,2])
                  .<= li.num_actions)
        @test all(time_series(MersenneTwister(seed), li, 3)
                  .<= li.num_actions)
        @test all(time_series(li, 3, init_actions) .<= li.num_actions)
        @test all(time_series(li, 3, [1,2]) .<= li.num_actions)
        @test @inferred(time_series(li, 3, init_actions, [1,2]) ==
                        [1 2 2; 1 1 2; 2 2 2])
        @test_throws ArgumentError series = time_series(li, 5, init_actions, [1, 1, 1])
    end

    @testset "Testing invalid local interaction instance" begin

        @testset "Non square adjacency matrix" begin
            non_square_adj = ones((2, 3))
            game = NormalFormGame(payoff_matrix)

            @test_throws ArgumentError li = LocalInteraction(game, non_square_adj)
            @test_throws ArgumentError li = LocalInteraction(payoff_matrix, non_square_adj)
        end

        @testset "Non square payoff matrix" begin
            non_square_pay = ones((2, 3))
            game = NormalFormGame((Player(ones((2, 3))), Player(ones((3, 2)))))

            @test_throws ArgumentError li = LocalInteraction(game, adj_matrix)
            @test_throws ArgumentError li = LocalInteraction(non_square_pay, adj_matrix)
        end
    end
end