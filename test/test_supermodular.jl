using Games

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

# Test script for supermodular

@testset "Testing supermodular.jl" begin
    
    @testset "is_supermodular for Player" begin 
        p1 = Player(zeros((2, 3)))
        coordination_game_matrix = [4 0; 3 2]
    	p2 = Player(coordination_game_matrix)
    	p3 = Player([0 1; 1 0])
    	@test is_supermodular(p1) == true
    	@test is_supermodular(p2) == true
    	@test is_supermodular(p3) == false
    end

    #@testset "is_supermodular for NormalFormGame" 
    #end
end
