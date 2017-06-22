# ------------------------- #
# Testing normal form games #
# ------------------------- #

# NOTE: We include `@inferred` at least once for each function name. We do
#       that multiple times for the same function if we have particular reason
#       to believe there might be a type stability with that function.

@testset "Testing game_theory/normal_form_game.jl" begin

    # Player #

    @testset "Player with 1 opponent" begin
        coordination_game_matrix = [4 0; 3 2]
        player = Player(coordination_game_matrix)

        @test @inferred(best_response(player, 2)) == 2
        @test @inferred(best_response(player, [1/2, 1/2])) == 2
        @test sort(@inferred(best_responses(player, [2/3, 1/3]))) ==
              sort([1, 2])
        @test best_response(
            player, [2/3, 1/3], tie_breaking="random"
            ) in [1, 2]
        @test @inferred(is_best_response(player, 1, 1))
        @test @inferred(is_best_response(player, [1/2, 1/2], [2/3, 1/3]))

        # Perturbed best response
        @test best_response(player, [2/3, 1/3], [0., 0.1]) == 2
        @test best_response(player, [2, 1], [0., 0.1]) == 2
    end

    @testset "Player with 2 opponents" begin
        payoffs_2opponents = Array{Int64}(2, 2, 2)
        payoffs_2opponents[:, 1, 1] = [3, 1]
        payoffs_2opponents[:, 1, 2] = [6, 0]
        payoffs_2opponents[:, 2, 1] = [4, 5]
        payoffs_2opponents[:, 2, 2] = [2, 7]
        player = @inferred Player(payoffs_2opponents)

        @test @inferred(payoff_vector(player, (1, 2))) == [6, 0]
        @test !(@inferred(is_best_response(player, 1, (2, 1))))
        @test @inferred(best_response(player, (2, 2))) == 2
        @test sort(@inferred(best_responses(player, ([3/7, 4/7], [1/2, 1/2])))) ==
              sort([1, 2])

        @test_throws MethodError best_response(player, (1, [1/2, 1/2]))
    end

    @testset "symmetric NormalFormGame with 2 players" begin
        coordination_game_matrix = [4 0; 3 2]
        g = @inferred(NormalFormGame(coordination_game_matrix))

        # NOTE: getindex(g, 1, 2) is equivalent to `g[1, 2]`. I use the former
        #       so we can test the type stability of the get index function
        @test @inferred(getindex(g, 1, 2)) == [0, 3]
        @test @inferred is_nash(g, (1, 1))
        @test @inferred is_nash(g, ([2/3, 1/3], [2/3, 1/3]))
    end

    # NormalFormGame #

    @testset "asymmetric NormalFormGame with 2 players" begin
        matching_pennies_bimatrix = Array{Float64}(2, 2, 2)
        matching_pennies_bimatrix[:, 1, 1] = [1, -1]
        matching_pennies_bimatrix[:, 1, 2] = [-1, 1]
        matching_pennies_bimatrix[:, 2, 1] = [-1, 1]
        matching_pennies_bimatrix[:, 2, 2] = [1, -1]
        g = @inferred(NormalFormGame(matching_pennies_bimatrix))

        @test g[2, 1] == [-1, 1]
        @test !(is_nash(g, (1, 1)))
        @test is_nash(g, ([1/2, 1/2], [1/2, 1/2]))
    end

    @testset "asymmetric NormalFormGame with 2 players" begin
        payoffs_2opponents = Array{Int64}(2, 2, 2)
        payoffs_2opponents[:, 1, 1] = [3, 1]
        payoffs_2opponents[:, 1, 2] = [6, 0]
        payoffs_2opponents[:, 2, 1] = [4, 5]
        payoffs_2opponents[:, 2, 2] = [2, 7]
        player = @inferred Player(payoffs_2opponents)
        g = @inferred NormalFormGame(tuple([player for i in 1:3]...))

        @test @inferred(getindex(g, 1, 1, 2)) == [6, 4, 1]
        @test @inferred is_nash(g, (1, 1, 1))
        @test @inferred !(is_nash(g, (1, 1, 2)))

        p = (1 + sqrt(65)) / 16
        @test is_nash(g, ([1-p, p], [1-p, p], [1-p, p]))
    end

    @testset "NormalFormGame input action sizes" begin
        g = @inferred NormalFormGame((2, 3, 4))

        @test @inferred(num_players(g)) == 3
        @test g.players[1].payoff_array == zeros((2, 3, 4))
        @test g.players[2].payoff_array == zeros((3, 4, 2))
        @test g.players[3].payoff_array == zeros((4, 2, 3))
    end

    @testset "NormalFormGame setindex" begin
        g = @inferred NormalFormGame((2, 2))
        g[1, 1] = [0, 10]
        g[1, 2] = [0, 10]
        g[2, 1] = [3, 5]
        g[2, 2] = [-2, 0]

        @test g.players[1].payoff_array == [0 0; 3 -2]
        @test g.players[2].payoff_array == [10 5; 10 0]
    end

    @testset "NormalFormGame constant payoffs" begin
        g = NormalFormGame((2, 2))

        @test @inferred is_nash(g, (1, 1))
        @test is_nash(g, (1, 2))
        @test is_nash(g, (2, 1))
        @test is_nash(g, (2, 2))
    end

    # Trivial cases with one player #

    @testset "Player with 0 opponents" begin
        payoffs = [0, 1]
        player = Player(payoffs)

        @test @inferred(payoff_vector(player, nothing)) == [0, 1]
        @test @inferred is_best_response(player, 2, nothing)
        @test @inferred(best_response(player, nothing)) == 2
    end

    @testset "NormalFormGame with 1 player" begin
        payoffs = [0, 1, 1]
        g = @inferred NormalFormGame(Player(payoffs))
        @test num_players(g) == 1
        @test g.players[1].payoff_array == [0, 1, 1]
        @test g[1] == 0
        @test is_nash(g, 2)
        @test !(is_nash(g, 1))
        @test is_nash(g, [0, 1/2, 1/2])

        g = NormalFormGame((2,))
        @test num_players(g) == 1
        @test g.players[1].payoff_array == zeros(2)
        g[1] = 10
        @test g.players[1].payoff_array == [10, 0]
    end

    # Invalid inputs #

    @testset "NormalFormGame invalid players shape inconsistent" begin
        p1 = Player(zeros((2, 3)))
        p2 = Player(zeros((2, 3)))
        @test_throws ArgumentError g = NormalFormGame((p1, p2))
    end

    @testset "NormalFormGame invalid players number inconsistent" begin
        p1 = Player(zeros((2, 2, 2)))
        p2 = Player(zeros((2, 2, 2)))
        @test_throws MethodError g = NormalFormGame((p1, p2))
    end

    @testset "NormalFormGame invalid nonsquare matrix" begin
        @test_throws ArgumentError g = NormalFormGame(zeros((2, 3)))
        @test_throws ArgumentError g = NormalFormGame(transpose([0 1 1]))
    end

    @testset "NormalFormGame invalid payoff profiles" begin
        @test_throws ArgumentError g = NormalFormGame(zeros((2, 2, 1)))
    end

    # Utility functions #

    @testset "pure2mixed" begin
        num_actions = 3
        pure_action = 1
        mixed_action = [1., 0., 0.]
        @test @inferred(pure2mixed(num_actions, pure_action)) == mixed_action
    end

end
