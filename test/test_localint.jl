# ------------------------- #
# Testing local interaction #
# ------------------------- #

# NOTE: We include `@inferred` at least once for each function name. We do
#       that multiple times for the same function if we have particular reason
#       to believe there might be a type stability with that function.


@testset "Testing localint.jl" begin
    payoff_matrix = [4 0; 2 3]
    adj_matrix = [0 1 3; 2 0 1; 3 2 0]
    init_actions = (1,1,2)

    @testset "Testing local interaction with simultaneous revision" begin
        li = LocalInteraction(payoff_matrix, adj_matrix)

        @test @inferred(play(li,init_actions) == (2,1,1))
        @test @inferred(time_series(li,3,init_actions) == [1 2 1; 1 1 2; 2 1 2])
    end

    @testset "Testing local interaction with sequencial revision" begin
        li = LocalInteraction(payoff_matrix, adj_matrix, AsynchronousRevision())

        @test @inferred(play(li,init_actions,1) == (2,1,2))
        @test @inferred(time_series(li,3,init_actions,[1,2]) ==
                        [1 2 2; 1 1 2; 2 2 2])
    end
end