
function state_EE(state::MixedDestabilizer, subsystem::Vector{Int})
    # Extract the stabilizer tableau and system parameters
    tab = state.tab
    n = tab.nqubits
    r = state.rank  # Number of independent stabilizers
    
    # Optimization: work with smaller subsystem for computational efficiency
    # Since we can use the complement for large subsystems
    if length(subsystem) <= div(n,2)
        AC = subsystem
    else
        AC = setdiff(1:n, subsystem)  # Use complement if subsystem > n/2
    end
    
    # Handle maximally mixed state (no stabilizers)
    # Entropy equals the subsystem size since all qubits are maximally mixed
    if r == 0
        return length(AC)
    end
    
    nsubsystem = length(AC)
    
    # Create binary matrix to represent how stabilizers act on subsystem AC
    # Rows: stabilizers, Columns: X and Z operators on each qubit in AC
    tab_part = falses(r, 2*nsubsystem)
    
    # Fill the matrix by checking each stabilizer's action on subsystem
    for i in 1:r
        for j in 1:nsubsystem
            # Get Pauli operator on qubit AC[j] from stabilizer i
            # Stabilizers are stored in rows (n+1) to (n+r) of the tableau
            op_ij = tab[i+n, AC[j]]
            
            # First half of columns: X components
            if op_ij[1]
                tab_part[i, j] = true
            end
            # Second half of columns: Z components
            if op_ij[2]
                tab_part[i, j + nsubsystem] = true
            end
        end
    end
    
    # Compute rank by row reduction and counting non-zero rows
    row_echelon!(tab_part)
    rnk = count(any(tab_part, dims=2))
    
    # Return entropy based on whether state is mixed or pure
    # Different formulas apply due to how information is encoded
    
	if r < length(AC)  # Mixed state
        return r - rnk # rnk - length(AC) + length(AC) - r = rnk - r
    else      # Pure state (r = n)
        return rnk - length(AC)
    end

end

function  row_echelon!(M::BitMatrix)
	K = size(M,1) # num generators (rows)
	N = div(size(M,2),2) # num qubits (cols)
	Ku = 1 # index of first row in the active region
	Nl = 1 # index of first column in the active region

	while ( Nl <= N && Ku <= K)
		num_x_ops = 0 # number of pauli x operators
		num_y_ops = 0 # number of pauli y operators
		num_z_ops = 0 # number of pauli z operators
		k1 = 0 # first row with pauli operator
		k1_op = "" # operator of k1
		k2 = 0 # first row with different pauli operator
		k2_op = "" # operator of k2

		for i in Ku:K
			if M[i,Nl] == 1 && M[i,Nl + N] == 0
				num_x_ops += 1
				if k1 == 0 && k1_op == ""
					k1 = i
					k1_op = "x"
				end

				if (k1_op == "y" || k1_op == "z") && k2 == 0 && k2_op == ""
					k2 = i
					k2_op = "x"
				end


			elseif M[i,Nl] == 0 && M[i,Nl + N] == 1
				num_z_ops += 1
				if k1 == 0 && k1_op == ""
					k1 = i
					k1_op = "z"
				end

				if (k1_op == "y" || k1_op == "x") && k2 == 0 && k2_op == ""
					k2 = i
					k2_op = "z"
				end

			elseif M[i,Nl] == 1 && M[i,Nl + N] == 1
				num_y_ops += 1
				if k1 == 0 && k1_op == ""
					k1 = i
					k1_op = "y"
				end

				if (k1_op == "x" || k1_op == "z") && k2 == 0 && k2_op == ""
					k2 = i
					k2_op = "y"
				end
			end
		end

		if (num_x_ops + num_y_ops + num_z_ops == 0)
			# No pauli operators in column Nl
			Nl += 1
		elseif (num_x_ops > 0 && (num_y_ops + num_z_ops == 0)) || (num_y_ops > 0 && (num_x_ops + num_z_ops == 0)) || (num_z_ops > 0 && (num_y_ops + num_x_ops == 0))
			# Only 1 kind of pauli operator
			one_kind_pauli(M,Ku,k1,Nl,K,N)
			Nl += 1
			Ku += 1
		else
			# There are at least 2 different kinds of pauli operators
			two_kind_pauli(M,Ku,k1,k2,Nl,K,N)
			Nl += 1
			Ku += 2
		end
	end


end

function one_kind_pauli(M::BitMatrix,Ku::Int,k::Int,Nl::Int,K::Int,N::Int)
	# Make row k top of active region by swapping (if necessary) row k with row Ku
	if k != Ku
		for j = 1:size(M,2)
			M[k,j],M[Ku,j] = M[Ku,j],M[k,j]
		end
	end

	# Multiply row Ku by all other rows that have same Pauli in column Nl
	for i in (Ku+1):K
		if M[Ku,Nl] == M[i,Nl] && M[Ku,Nl+N] == M[i,Nl+N]
			for cn = 1:N
				M[i,cn] = (M[i,cn] + M[Ku,cn])%2
				M[i,cn+N] = (M[i,cn+N] + M[Ku,cn+N])%2
			end
		end
	end
end

function two_kind_pauli(M::BitMatrix,Ku::Int,k1::Int,k2::Int,Nl::Int,K::Int,N::Int)
	# Make row k1 top of active region by swapping (if necessary) row k1 with row Ku
	if k1 != Ku
		for j = 1:size(M,2)
			M[k1,j],M[Ku,j] = M[Ku,j],M[k1,j]
		end
	end

	# Make row k2 second top of active region by swapping (if necessary) row k2 with row Ku+1
	if k2 != Ku+1
		for j = 1:size(M,2)
			M[k2,j],M[Ku+1,j] = M[Ku+1,j],M[k2,j]
		end
	end

	for i in (Ku+2):K
		if M[i,Nl] == 0 && M[i,Nl+N] == 0
			# Do nothing
		elseif M[Ku,Nl] == M[i,Nl] && M[Ku,Nl+N] == M[i,Nl+N]
			for cn = 1:N
				M[i,cn] = (M[i,cn] + M[Ku,cn])%2
				M[i,cn+N] = (M[i,cn+N] + M[Ku,cn+N])%2
			end
		elseif M[Ku+1,Nl] == M[i,Nl] && M[Ku+1,Nl+N] == M[i,Nl+N]
			for cn = 1:N
				M[i,cn] = (M[i,cn] + M[Ku+1,cn])%2
				M[i,cn+N] = (M[i,cn+N] + M[Ku+1,cn+N])%2
			end
		else
			for cn = 1:N
				M[i,cn] = (M[i,cn] + M[Ku,cn])%2
				M[i,cn+N] = (M[i,cn+N] + M[Ku,cn+N])%2
			end
			for cn = 1:N
				M[i,cn] = (M[i,cn] + M[Ku+1,cn])%2
				M[i,cn+N] = (M[i,cn+N] + M[Ku+1,cn+N])%2
			end
		end
	end
end