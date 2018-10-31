# ------------------------- #
# Testing normal form games #
# ------------------------- #

# NOTE: We include `@inferred` at least once for each function name. We do
#       that multiple times for the same function if we have particular reason
#       to believe there might be a type stability with that function.

@testset "Testing normal_form_game.jl" begin

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
        @test_throws ArgumentError best_response(
            player, [2/3, 1/3], tie_breaking="0"
            )
        @test @inferred(is_best_response(player, 1, 1))
        @test @inferred(is_best_response(player, [1/2, 1/2], [2/3, 1/3]))

        # Perturbed best response
        @test best_response(player, [2/3, 1/3], [0., 0.1]) == 2
        @test best_response(player, [2, 1], [0., 0.1]) == 2
    end

    @testset "Player with 2 opponents" begin
        payoffs_2opponents = Array{Int64}(undef, 2, 2, 2)
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

    @testset "repr(Player)" begin
        A = [1 2; 3 4]
        player = Player(A)
        r = repr("text/plain", A)
        @test repr(player) ==
            replace(r, string(typeof(A)) =>
                       split(string(typeof(player)), ".")[end])
    end

    # NormalFormGame #

    @testset "symmetric NormalFormGame with 2 players" begin
        coordination_game_matrix = [4 0; 3 2]
        g = @inferred(NormalFormGame(coordination_game_matrix))

        # NOTE: getindex(g, 1, 2) is equivalent to `g[1, 2]`. I use the former
        #       so we can test the type stability of the get index function
        @test @inferred(getindex(g, 1, 2)) == [0, 3]
        @test @inferred(getindex(g, CartesianIndex(1, 2))) == [0, 3]
        @test @inferred is_nash(g, (1, 1))
        @test @inferred is_nash(g, ([2/3, 1/3], [2/3, 1/3]))
    end

    @testset "asymmetric NormalFormGame with 2 players" begin
        matching_pennies_bimatrix = Array{Float64}(undef, 2, 2, 2)
        matching_pennies_bimatrix[:, 1, 1] = [1, -1]
        matching_pennies_bimatrix[:, 1, 2] = [-1, 1]
        matching_pennies_bimatrix[:, 2, 1] = [-1, 1]
        matching_pennies_bimatrix[:, 2, 2] = [1, -1]
        g = @inferred(NormalFormGame(matching_pennies_bimatrix))

        @test g[2, 1] == [-1, 1]
        @test !(is_nash(g, (1, 1)))
        @test is_nash(g, ([1/2, 1/2], [1/2, 1/2]))
    end

    @testset "asymmetric NormalFormGame with 3 players" begin
        payoffs_2opponents = Array{Int64}(undef, 2, 2, 2)
        payoffs_2opponents[:, 1, 1] = [3, 1]
        payoffs_2opponents[:, 1, 2] = [6, 0]
        payoffs_2opponents[:, 2, 1] = [4, 5]
        payoffs_2opponents[:, 2, 2] = [2, 7]
        player = @inferred Player(payoffs_2opponents)
        g = @inferred NormalFormGame(tuple([player for i in 1:3]...))

        @test @inferred(getindex(g, 1, 1, 2)) == [6, 4, 1]
        @test @inferred(getindex(g, CartesianIndex(1, 1, 2))) == [6, 4, 1]
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
        g[CartesianIndex(2, 2)] = [-2, 0]

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
        @test is_dominated(player, 1)
        @test !is_dominated(player, 2)
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

        A = [0, 1, 1]
        A = reshape(A, (size(A)..., 1))
        @test_throws ArgumentError g = NormalFormGame(A)
    end

    @testset "NormalFormGame invalid payoff profiles" begin
        @test_throws ArgumentError g = NormalFormGame(zeros((2, 2, 1)))
    end

    @testset "NormalFormGame empty tuple" begin
        @test_throws ArgumentError g = NormalFormGame(tuple())
    end

    @testset "payoff_vector empty tuple" begin
        p1 = Player(zeros((2, 2, 2)))
        @test_throws ArgumentError payoff_vector(p1, tuple())
    end

    # Utility functions #

    @testset "pure2mixed" begin
        num_actions = 3
        pure_action = 1
        mixed_action = [1., 0., 0.]
        @test @inferred(pure2mixed(num_actions, pure_action)) == mixed_action
    end

    # Pareto efficiency & Pareto dominance #

    @testset "Tests on Pareto efficiency and dominance" begin
        coordination_game_matrix = [4 0;
                                    3 2]

        equal_po_p1_bimatrix = Array{Float64}(undef, 2, 2, 2)
        equal_po_p1_bimatrix[1, 1, :] = [1, -1]
        equal_po_p1_bimatrix[1, 2, :] = [1, 1]
        equal_po_p1_bimatrix[2, 1, :] = [1, 1]
        equal_po_p1_bimatrix[2, 2, :] = [1, -1]

        three_p_equal_po_array = Array{Int64}(undef, 2, 2, 2)
        three_p_equal_po_array[:, :, 1] = [2 0; 0 2]
        three_p_equal_po_array[:, :, 2] = [2 0; 0 2]

        p1 = p2 = p3 = Player(three_p_equal_po_array)

        games_dict = [NormalFormGame(coordination_game_matrix),
                      NormalFormGame(equal_po_p1_bimatrix),
                      NormalFormGame((p1, p2, p3))]

        act_prof_dict = [[(1, 1), (1, 2), (2, 1), (2, 2)],
                         [(1, 1), (1, 2), (2, 1), (2, 2)],
                         [(1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1),
                          (1, 1, 2), (2, 1, 2), (1, 2, 2), (2, 2, 2)]]

        @testset "Testing is_pareto_efficient" begin
            output_dict = [[true, false, false, false],
                           [false, true, true, false],
                           [true, false, false, false, false, false, false,
                            true]]

            for i = 1:length(games_dict)
                for j in 1:length(act_prof_dict[i])
                    @test @inferred is_pareto_efficient(games_dict[i],
                                                       act_prof_dict[i][j]) ==
                                    output_dict[i][j]
                end
            end
        end

        @testset "Testing Pareto dominance" begin
            output_dict = [[true, false, false, false],
                           [false, false, false, false],
                           [false, false, false, false, false, false, false,
                            false]]

            for i = 1:length(games_dict)
                for j in 1:length(act_prof_dict[i])
                    @test @inferred is_pareto_dominant(games_dict[i],
                                                       act_prof_dict[i][j]) ==
                                    output_dict[i][j]
                end
            end
        end

        @testset "Test is_dominated" begin
            coordination_game_matrix = [4 0; 3 2]
            player = Player(coordination_game_matrix)
            for action = 1:num_actions(player)
                @test !is_dominated(player, action)
            end

            payoffs_2opponents = Array{Int64}(undef, 2, 2, 2)
            payoffs_2opponents[:, 1, 1] = [3, 1]
            payoffs_2opponents[:, 1, 2] = [6, 0]
            payoffs_2opponents[:, 2, 1] = [4, 5]
            payoffs_2opponents[:, 2, 2] = [2, 7]
            player = Player(payoffs_2opponents)

            for i = 1:num_actions(player)
                @test !is_dominated(player, i)
            end

        end

        @testset "Test player corner cases" begin
            n, m = 3, 4
            player = Player(zeros((n, m)))
            for action = 1:n
                @test is_best_response(player, action, ones(m) * 1/m)
                @test !is_dominated(player, action)
            end

            e = 1e-8
            player = Player([-e -e;
                             1 -1;
                             -1 1])
            action = 1
            @test is_best_response(player, action, [1/2, 1/2], tol=e)
            @test !is_best_response(player, action, [1/2, 1/2], tol=e/2)
            @test !is_dominated(player, action, tol=e+1e-16)
            @test is_dominated(player, action, tol=e/2)
        end
    end

end
