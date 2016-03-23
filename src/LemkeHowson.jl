#=
Katrina Evtimova, kve216@nyu.edu
Code adapted from: http://www.mathworks.com/matlabcentral/fileexchange/44279-lemke-howson-algorithm-for-2-player-games/content/LemkeHowson.m
=#

#=
function LemkeHowson(game::NormalFormGame, k0::Int, maxPivots::Int)
# Input: two player game
# Output: Nash Equilibrium - array with two strategies for player 1 and 2, respectively

	# Re-scale payoffs to positive

	# Tableaus

	# Row Labels

	# Complimentary pivoting

	# Pivoting Loop

	# Calculating Nash

end
=#

function pivot(player::Player{2}, r::Int, s::Int)
	m, _ = size(player.payoff_array)
	M = copy(player.payoff_array)

	for i = 1:m
		if i != r
			M[i, :] = player.payoff_array[i,:] - player.payoff_array[i,s]/player.payoff_array[r,s]*player.payoff_array[r,:]
		else
			continue
		end
	end

	return M
end