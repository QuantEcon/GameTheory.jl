#=
Katrina Evtimova, kve216@nyu.edu
Code adapted from: http://www.mathworks.com/matlabcentral/fileexchange/44279-lemke-howson-algorithm-for-2-player-games/content/LemkeHowson.m
=#


function LemkeHowson(game::NormalFormGame, k0::Int=1, max_pivots::Int=500000)
# Input: two player game
# Output: Nash Equilibrium - array with two strategies for player 1 and 2, respectively

	# Check input 
	if length(g.players) != 2
		throw(ArgumentError("2-player game expected"))
    end

    # Get payoff arrays
    A = game.players[1].payoff_array
    B = game.players[2].payoff_array
    m,n = size(A)

	# Re-scale payoffs to positive
	min_val = min(minimum(A),minimum(B))
	if min_val <=0 
		A = A - min_val + 1
		B = B - min_val + 1
	end

	# Tableaus
	tableaus = []
	tab_1 = [B' eye(n) ones(n,1)]
	tab_2 = [eye(m) A ones(m,1)]
	push!(tableaus, tab_1, tab_2)

	# Row Labels
	row_labels = []
	row_lab_1 = [m+1:m+n]
	row_lab_2 = [1:m]
	push!(row_labels, row_lab_1, row_lab_2)

	# Complimentary pivoting
	k = k0
	if m <= n
		player = 1
	else 
		player = 2
	end

	# Pivoting Loop
	num_piv = 0
	while num_piv < max_piv
		num_piv += 1

		# Pick tableau
		tableau = tableaus[player]
		m_, _ = size(tableau)

		# Find pivot row
		max_ = 0 
		ind = -1

		for i = 1:m
			ratio = tableau[i,k]/tableau[i, m+n+1]
			if ratio > max_
				ind = i
				max_ = ratio
			end
		end

		if max_ > 0
			tableaus[player] = pivot(tableau, ind, k)
		else
			break
		end

		# swap labels
		row_labels[player][ind], k = k, row_labels[player][ind]

		# Break if entering variable equals k0
		if k == k0
			break
		end

		# Alternate players
		player = 3 - player

		# Return error if max_pivot iterations are reached
		if num_piv == max_pivots
			throw(ErrorException("Max number of iterations reached."))
		end

	end

	# Calculating Nash
	NashEquil = []

	for pl = [1:2]
		x = zeros(size(A)[pl], 1)
		rows = row_labels[pl]
		tableau = tableaus[pl]

		for i = [1:length(rows)]
			if pl == 1 && rows[i] == size(A)[1]
				x[rows[i]] = tableau[i,m+n+1]/tableau[i,rows[i]]
			elseif pl == 2 && rows[i] > size(A)[1]
				x[rows[i] - size(A)[1]] = tableau[i,m+n-1]/tableau[i,rows[i]]
			end
		end

		push!(NashEquil, x/sum(x))
	end

	return NashEquil

end


function pivot(tableau, r::Int, s::Int)
	m, _ = size(tableau)
	M = copy(tableau)

	for i = 1:m
		if i != r
			M[i, :] = tableau[i,:] - tableau[i,s]/tableau[r,s]*tableau[r,:]
		else
			continue
		end
	end

	return M
end