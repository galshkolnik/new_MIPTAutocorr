"""
Apply a single qubit Clifford operation to a quantum state.

Parameters:
- clifford_index: Integer between 1 and 24 representing the Clifford operation
- qubit_index: Index of the qubit to apply the operation to
- state: QuantumClifford stabilizer state to modify

Returns:
- Modified state after applying the Clifford operation
"""
function apply_single_qubit_clifford!(clifford_index::Int, qubit_index::Int, state)
    if clifford_index < 1 || clifford_index > 24
        throw(ArgumentError("clifford_index must be between 1 and 24"))
    end
    
    # Define the 24 single qubit Cliffords
    # Each entry is a list of gates to apply in sequence
    clifford_operations = [
        # Paulis
        [],                          # 1: I (Identity)
        [sX(qubit_index)],          # 2: X
        [sY(qubit_index)],          # 3: Y
        [sY(qubit_index), sX(qubit_index)],  # 4: Y, X
        
        # 2π/3 rotations
        [sSQRTX(qubit_index), sSQRTY(qubit_index)],                # 5: X/2, Y/2
        [sSQRTX(qubit_index), sInvSQRTY(qubit_index)],             # 6: X/2, -Y/2
        [sInvSQRTX(qubit_index), sSQRTY(qubit_index)],             # 7: -X/2, Y/2
        [sInvSQRTX(qubit_index), sInvSQRTY(qubit_index)],          # 8: -X/2, -Y/2
        [sSQRTY(qubit_index), sSQRTX(qubit_index)],                # 9: Y/2, X/2
        [sSQRTY(qubit_index), sInvSQRTX(qubit_index)],             # 10: Y/2, -X/2
        [sInvSQRTY(qubit_index), sSQRTX(qubit_index)],             # 11: -Y/2, X/2
        [sInvSQRTY(qubit_index), sInvSQRTX(qubit_index)],          # 12: -Y/2, -X/2
        
        # π/2 rotations
        [sSQRTX(qubit_index)],         # 13: X/2
        [sInvSQRTX(qubit_index)],      # 14: -X/2
        [sSQRTY(qubit_index)],         # 15: Y/2
        [sInvSQRTY(qubit_index)],      # 16: -Y/2
        [sInvSQRTX(qubit_index), sSQRTY(qubit_index), sSQRTX(qubit_index)],  # 17: -X/2, Y/2, X/2
        [sInvSQRTX(qubit_index), sInvSQRTY(qubit_index), sSQRTX(qubit_index)],  # 18: -X/2, -Y/2, X/2
        
        # Hadamard-like
        [sX(qubit_index), sSQRTY(qubit_index)],                    # 19: X, Y/2
        [sX(qubit_index), sInvSQRTY(qubit_index)],                 # 20: X, -Y/2
        [sY(qubit_index), sSQRTX(qubit_index)],                    # 21: Y, X/2
        [sY(qubit_index), sInvSQRTX(qubit_index)],                 # 22: Y, -X/2
        [sSQRTX(qubit_index), sSQRTY(qubit_index), sSQRTX(qubit_index)],  # 23: X/2, Y/2, X/2
        [sInvSQRTX(qubit_index), sSQRTY(qubit_index), sInvSQRTX(qubit_index)]   # 24: -X/2, Y/2, -X/2
    ]
    
    # Apply the sequence of gates for the chosen Clifford operation
    operations = clifford_operations[clifford_index]
    for op in operations
        apply!(state, op)
    end
    
    return state
end

"""
Create a lookup table for the 24 single qubit Cliffords for reference.
"""
function get_clifford_description(clifford_index::Int)
    descriptions = [
        "I",                    # 1
        "X",                    # 2
        "Y",                    # 3
        "Y, X",                 # 4
        "X/2, Y/2",            # 5
        "X/2, -Y/2",           # 6
        "-X/2, Y/2",           # 7
        "-X/2, -Y/2",          # 8
        "Y/2, X/2",            # 9
        "Y/2, -X/2",           # 10
        "-Y/2, X/2",           # 11
        "-Y/2, -X/2",          # 12
        "X/2",                 # 13
        "-X/2",                # 14
        "Y/2",                 # 15
        "-Y/2",                # 16
        "-X/2, Y/2, X/2",      # 17
        "-X/2, -Y/2, X/2",     # 18
        "X, Y/2",              # 19
        "X, -Y/2",             # 20
        "Y, X/2",              # 21
        "Y, -X/2",             # 22
        "X/2, Y/2, X/2",       # 23
        "-X/2, Y/2, -X/2"      # 24
    ]
    
    if clifford_index < 1 || clifford_index > 24
        throw(ArgumentError("clifford_index must be between 1 and 24"))
    end
    
    return descriptions[clifford_index]
end


"""
Apply a two-qubit Clifford operation to a quantum state.

Parameters:
- qubit1: Index of the first qubit
- qubit2: Index of the second qubit 
- gate_index: Integer between 0 and 11519 representing the two-qubit Clifford operation
- state: QuantumClifford stabilizer state to modify

Returns:
- Modified state after applying the two-qubit Clifford operation
"""
function apply_random_two_qubit_clifford!(qubit1::Int, qubit2::Int, state; gate_index::Int=-1)
    # Randomly choose gate index between 0 and 11519
    if gate_index==-1; gate_index = rand(0:11519); end

    if gate_index < 0 || gate_index > 11519
        throw(ArgumentError("gate must be between 0 and 11519"))
    end

    ind1 = mod(gate_index, 576)
    # Apply initial single-qubit gates
    first_gate_index = div(ind1, 24) + 1
    second_gate_index = mod(ind1, 24) + 1
    apply_single_qubit_clifford!(first_gate_index, qubit1, state)
    apply_single_qubit_clifford!(second_gate_index, qubit2, state)
    
    if gate_index < 576
        # Case 1: Identity-class, single qubit gates only
        
    elseif gate_index < 5760
        # Case 2: CNOT-like class
        
        # Apply CNOT gate
        apply!(state, sCNOT(qubit1, qubit2))
        
        # Apply final single-qubit gates from S1 sets
        ind2 = div(gate_index - 576, 576)
        i1 = div(ind2, 3) + 1
        i2 = mod(ind2, 3) + 1
        
        # Apply final gate on first qubit
        if i1 == 2
            apply_single_qubit_clifford!(8, qubit1, state)  # "-X/2, -Y/2"
        elseif i1 == 3
            apply_single_qubit_clifford!(9, qubit2, state)  # "Y/2, X/2"
        end
        # i1 == 1: do nothing (identity)
        
        # Apply final gate on second qubit
        if i2 == 2
            apply_single_qubit_clifford!(8, qubit2, state)  # "-X/2, -Y/2"
        elseif i2 == 3
            apply_single_qubit_clifford!(9, qubit2, state)  # "Y/2, X/2"
        end
        # i2 == 1: do nothing (identity)
        
    elseif gate_index < 10944
        # Case 3: iSWAP-like class
        gate_index_adjusted = gate_index - 5760
        
        # Apply iSWAP gate
        apply!(state, sISWAP(qubit1, qubit2))
        
        # Apply final single-qubit gates from S1 sets
        ind2 = div(gate_index_adjusted, 576)
        i1 = div(ind2, 3) + 1
        i2 = mod(ind2, 3) + 1
        
        # Apply final gate on first qubit
        if i1 == 2
            apply_single_qubit_clifford!(8, qubit1, state)  # "-X/2, -Y/2"
        elseif i1 == 3
            apply_single_qubit_clifford!(9, qubit1, state)  # "Y/2, X/2"
        end
        
        # Apply final gate on second qubit
        if i2 == 2
            apply_single_qubit_clifford!(8, qubit2, state)  # "-X/2, -Y/2"
        elseif i2 == 3
            apply_single_qubit_clifford!(9, qubit2, state)  # "Y/2, X/2"
        end
        
    else  # 10944 <= gate_index < 11520
        # Case 4: SWAP-like class
        
        # Apply SWAP gate
        apply!(state, sSWAP(qubit1, qubit2))
    end
    
    return state
end

"""
Apply a two-qubit dual-unitary operation to a quantum state.

Parameters:
- qubit1: Index of the first qubit
- qubit2: Index of the second qubit 
- gate_index: Integer between 0 and 5759 representing the two-qubit dual-unitary operation
- state: QuantumClifford stabilizer state to modify

Returns:
- Modified state after applying the two-qubit dual unitary operation
"""
function apply_random_two_qubit_dual_unitary!(qubit1::Int, qubit2::Int, state; gate_index::Int=-1)
    # Randomly choose gate index between 0 and 5759
    if gate_index==-1; gate_index = rand(0:5759); end

    if gate_index < 0 || gate_index > 5759
        throw(ArgumentError("gate must be between 0 and 5759"))
    end

    ind1 = mod(gate_index, 576)
    # Apply initial single-qubit gates
    first_gate_index = div(ind1, 24) + 1
    second_gate_index = mod(ind1, 24) + 1
    apply_single_qubit_clifford!(first_gate_index, qubit1, state)
    apply_single_qubit_clifford!(second_gate_index, qubit2, state)
    
    if gate_index < 5184
        # Case 1: iSWAP-like class
        
        # Apply iSWAP gate
        apply!(state, sISWAP(qubit1, qubit2))
        
        # Apply final single-qubit gates from S1 sets
        ind2 = div(gate_index, 576)
        i1 = div(ind2, 3) + 1
        i2 = mod(ind2, 3) + 1
        
        # Apply final gate on first qubit
        if i1 == 2
            apply_single_qubit_clifford!(8, qubit1, state)  # "-X/2, -Y/2"
        elseif i1 == 3
            apply_single_qubit_clifford!(9, qubit1, state)  # "Y/2, X/2"
        end
        
        # Apply final gate on second qubit
        if i2 == 2
            apply_single_qubit_clifford!(8, qubit2, state)  # "-X/2, -Y/2"
        elseif i2 == 3
            apply_single_qubit_clifford!(9, qubit2, state)  # "Y/2, X/2"
        end
        
    else  # 5184 <= gate_index < 5760
        # Case 2: SWAP-like class
        
        # Apply SWAP gate
        apply!(state, sSWAP(qubit1, qubit2))
    end
    
    return state
end

"""
Perform a Z-basis measurement on a specified qubit.

Parameters:
- qubit_index: Index of the qubit to measure
- state: QuantumClifford stabilizer state to modify

Returns:
- state: Modified state after the measurement (projected onto the measurement outcome)
- outcome: Measurement outcome (+1 for |0⟩, -1 for |1⟩)
"""
function measure_z!(qubit_index::Int, state; meas_outcome::Int=0)
    if state.rank>0 # If the state is not maximally mixed, we can measure directly
        # Perform the measurement and get the outcome
        newstate, anticomindex, meas_res = projectZ!(state, qubit_index)
        is_det = true
        if isnothing(meas_res)
            meas_res = meas_outcome==0 ? rand((-1,1)) : meas_outcome
            m = meas_res==1 ? 0x0 : 0x2
            phases(newstate)[anticomindex] = m
            is_det = false
        end
    elseif state.rank==0 # If the state is maximally mixed, we need to use project!
        L=state.tab.nqubits
        Zinds=falses(L); Zinds[qubit_index] = true
        phase = meas_outcome==0 ? rand((0x0,0x2)) : meas_outcome==1 ? 0x0 : 0x2
        stab_q= Stabilizer([phase],transpose(zeros(Bool,L)), transpose(Zinds))  # Create a stabilizer with Z on the qubit
        newstate=MixedDestabilizer(stab_q)  # Create a new state with this stabilizer
        meas_res = phase==0x0 ? 1 : -1
        is_det = false
    else
        throw(ArgumentError("Invalid state for measurement"))
    end
    
    return newstate, meas_res, is_det # The measurement outcome
end