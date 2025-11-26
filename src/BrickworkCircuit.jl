"""
Create initial state for the monitored quantum circuit.

Parameters:
- L: Number of qubits
- state_type: Symbol indicating the type of initial state
  - :product_0: Product state |0⟩^L
  - :product_rand: Product state with random Z-basis states |0⟩ or |1⟩ 
  - :Bell_pairs: Bell pairs (|00⟩ + |11⟩) for each pair of qubits
  - :Bell_pairs_nnn: Bell pairs (|00⟩ + |11⟩) for next-nearest neighboring qubits
  - :maximally_entangled: Maximally entangled state
  - :fully_mixed: Fully mixed state (maximum entropy)

Returns:
- MixedDestabilizer state ready for circuit simulation
"""
function create_initial_state(L::Int, state_type::Symbol)
    if state_type == :product_0
        # All qubits in |0⟩
        return MixedDestabilizer(one(Stabilizer, L))
        
    elseif state_type == :product_rand
        # Start with |0⟩^L and randomly flip some qubits to |1⟩
        state = MixedDestabilizer(one(Stabilizer, L))
        for i in 1:L
            if rand() < 0.5  # 50% chance to flip each qubit
                apply!(state, sX(i))  # Apply X gate: |0⟩ → |1⟩
            end
        end
        return state

    elseif state_type == :Bell_pairs
        # Create Bell pairs: |00⟩ + |11⟩ for each pair of qubits
        if L % 2 != 0
            throw(ArgumentError("L must be even for Bell pairs"))
        end
        state = MixedDestabilizer(one(Stabilizer, L))
        for i in 1:2:L-1
            apply!(state, sHadamard(i))  # Apply Hadamard to first qubit of the pair
            apply!(state, sCNOT(i, i+1))  # Apply CNOT to create Bell pair
        end
        return state
    elseif state_type == :Bell_pairs_nnn
        # Create Bell pairs: |00⟩ + |11⟩ for each next-nearest neighbors (1,3),(2,4),(5,7),...
        if L % 2 != 0
            throw(ArgumentError("L must be even for Bell pairs"))
        end
        state = MixedDestabilizer(one(Stabilizer, L))
        for i in 1:2:L-1
            if mod1(i,4)>2
                apply!(state, sHadamard(i-2))  # Apply Hadamard to first qubit of the pair
                apply!(state, sCNOT(i-2, i))  # Apply CNOT to create Bell pair
            end
        end
        return state
        
    elseif state_type == :fully_mixed
        # Fully mixed state - rank 0 mixed state (no stabilizer generators)
        base_stab = zero(Stabilizer, L)  # Create L-qubit tableau
        return MixedDestabilizer(tab(base_stab), 0)  # Set rank to 0 for maximum mixing
    elseif state_type == :maximally_entangled
        # Maximally entangled state
        zs_op = zeros(Bool, L, L); xs_op= zeros(Bool, L, L)
        for i in 1:L
            zs_op[i, L+1-i] = true  # Z on qubit i
            xs_op[i, i] = true  # X on qubit i
        end
        phases=rand((0x0,0x2),L)
        stab_q= Stabilizer(phases,xs_op,zs_op) 
        return MixedDestabilizer(stab_q)
        
    else
        throw(ArgumentError("Invalid state_type. Must be :product_0, :product_rand, :Bell_pairs, :Bell_pairs_nnn, :fully_mixed or :maximally_entangled"))
    end
end


"""
Simulate a brickwork monitored quantum circuit.

Parameters:
- L: Number of qubits
- T: Number of time layers  
- p: Measurement probability (0 ≤ p ≤ 1)
- initial_state_type: Type of initial state (:product_0, :product_rand, :fully_mixed)
- is_pbc: Boolean for periodic boundary conditions

Returns:
- final_state: The quantum state after T layers
- measurement_record: Array of measurement outcomes for each layer
"""
function brickwork_circuit(L::Int, T::Int, p::Float64, state, is_pbc::Bool; extract_gates::Bool=false, measure_first_qubit::Bool=true, unitaries_type::Symbol=:Cliffords)
    # Validate inputs
    if !(0 ≤ p ≤ 1)
        throw(ArgumentError("Measurement probability p must be between 0 and 1"))
    end
    if !(unitaries_type in [:Cliffords, :DualUnitaries])
        throw(ArgumentError("Invalid unitaries_type. Must be :Cliffords or :DualUnitaries"))
    end
    
    # Pre-compute bond lists (create once, use many times)
    # === ODD BONDS: [1,2], [3,4], [5,6], ... ===
    odd_bonds = []
    for i in 1:div(L, 2)
        push!(odd_bonds, (2*i-1, 2*i))
    end
    # Add periodic boundary condition for odd bonds if L is odd
    if is_pbc && (L % 2 == 1)
        push!(odd_bonds, (L, 1))
    end
    
    # === EVEN BONDS: [2,3], [4,5], [6,7], ... ===
    even_bonds = []
    for i in 1:div(L-1, 2)
        push!(even_bonds, (2*i, 2*i+1))
    end
    # Add periodic boundary condition for even bonds if L is even
    if is_pbc && (L % 2 == 0)
        push!(even_bonds, (L, 1))
    end
    
    # Track measurement outcomes: [time_step][qubit] -> outcome (or nothing if not measured)
    measurement_record = []
    gates_record = []
    
    # Main circuit loop
    for t in 1:T
        # === EVEN BONDS LAYER ===
        if extract_gates
            gates_realization_even = Vector{Union{Nothing, Int}}(nothing, length(even_bonds))
        end
        # Apply random Clifford gates on even bonds
        for (bond_ind,(q1, q2)) in enumerate(even_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            
            if extract_gates; gates_realization_even[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_even); end
        
        # === FIRST MEASUREMENT LAYER ===
        measurement_outcomes_1 = Vector{Union{Nothing, Int}}(nothing, L)
        for qubit in 1:L
            if qubit == 1 && (!measure_first_qubit)
                continue  # Skip measurement of first qubit if not required
            elseif rand() < p  # Measure with probability p
                state, outcome, is_det = measure_z!(qubit, state)
                measurement_outcomes_1[qubit] = outcome
            end
        end
        push!(measurement_record, measurement_outcomes_1)

        # === ODD BONDS LAYER ===
        if extract_gates
            gates_realization_odd = Vector{Union{Nothing, Int}}(nothing, length(odd_bonds))
        end
        # Apply random Clifford gates on odd bonds
        for (bond_ind,(q1, q2)) in enumerate(odd_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            if extract_gates; gates_realization_odd[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_odd); end
        
        # === SECOND MEASUREMENT LAYER ===
        measurement_outcomes_2 = Vector{Union{Nothing, Int}}(nothing, L)
        if t<T
            for qubit in 1:L
                if qubit == 1 && (!measure_first_qubit)
                    continue  # Skip measurement of first qubit if not required
                elseif rand() < p  # Measure with probability p
                    state, outcome, is_det = measure_z!(qubit, state)
                    measurement_outcomes_2[qubit] = outcome
                end
            end
        end
        push!(measurement_record, measurement_outcomes_2)
    end
    
    if extract_gates
        return state, measurement_record, gates_record
    else
        return state, measurement_record
    end
end


"""
Simulate a brickwork monitored quantum circuit with given gates and measurements.
"""
function brickwork_circuit(L::Int, T::Int, state, is_pbc::Bool; measurement_gates::Vector{Any}, unitary_gates::Vector{Any}, unitaries_type::Symbol=:Cliffords)

    # Validate inputs
    if !(length(unitary_gates) == 2*T && length(measurement_gates) == 2*T)
        throw(ArgumentError("Invalid gate configuration: must have 2T layers of gates"))
    end
        if !(unitaries_type in [:Cliffords, :DualUnitaries])
        throw(ArgumentError("Invalid unitaries_type. Must be :Cliffords or :DualUnitaries"))
    end
    
    # Pre-compute bond lists (create once, use many times)
    # === ODD BONDS: [1,2], [3,4], [5,6], ... ===
    odd_bonds = []
    for i in 1:div(L, 2)
        push!(odd_bonds, (2*i-1, 2*i))
    end
    # Add periodic boundary condition for odd bonds if L is odd
    if is_pbc && (L % 2 == 1)
        push!(odd_bonds, (L, 1))
    end
    
    # === EVEN BONDS: [2,3], [4,5], [6,7], ... ===
    even_bonds = []
    for i in 1:div(L-1, 2)
        push!(even_bonds, (2*i, 2*i+1))
    end
    # Add periodic boundary condition for even bonds if L is even
    if is_pbc && (L % 2 == 0)
        push!(even_bonds, (L, 1))
    end
    
    # Track measurement outcomes: [time_step][qubit] -> outcome (or nothing if not measured)
    measurement_record = []
    
    # Main circuit loop
    for t in 1:T
        # === EVEN BONDS LAYER ===
        # Apply random Clifford gates on even bonds
        for (jj,(q1, q2)) in enumerate(even_bonds)
            gate_index = unitary_gates[2*t-1][jj]
            if unitaries_type == :Cliffords
                # Use Clifford gates
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
        end
        
        # === FIRST MEASUREMENT LAYER ===
        measurement_outcomes_1 = Vector{Union{Nothing, Int}}(nothing, L)
        for qubit in 1:L
            if !isnothing(measurement_gates[2*t-1][qubit]) # Measure this qubit
                state, outcome, is_det = measure_z!(qubit, state, meas_outcome=measurement_gates[2*t-1][qubit])
                measurement_outcomes_1[qubit] = outcome
            end
        end
        push!(measurement_record, measurement_outcomes_1)
        
        # === ODD BONDS LAYER ===
        # Apply random Clifford gates on odd bonds
        for (jj,(q1, q2)) in enumerate(odd_bonds)
            gate_index = unitary_gates[2*t][jj]
            if unitaries_type == :Cliffords
                # Use Clifford gates
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
        end
        
        # === SECOND MEASUREMENT LAYER ===
        measurement_outcomes_2 = Vector{Union{Nothing, Int}}(nothing, L)
        if t<T
            for qubit in 1:L
                if !isnothing(measurement_gates[2*t][qubit]) # Measure this qubit
                    state, outcome, is_det = measure_z!(qubit, state, meas_outcome=measurement_gates[2*t][qubit])
                    measurement_outcomes_2[qubit] = outcome
                end
            end
        end
        push!(measurement_record, measurement_outcomes_2)
    end
    
    return state, measurement_record
end

"""
Allows for quenched disorder.
"""
function brickwork_circuit(L::Int, T::Int, p_ar::Vector{Float64}, state, is_pbc::Bool; extract_gates::Bool=false, measure_first_qubit::Bool=true, unitaries_type::Symbol=:Cliffords)
    # Validate inputs
    if !(length(p_ar)==L)
        throw(ArgumentError("probability vector p_ar must have length L"))
    end
    if !(unitaries_type in [:Cliffords, :DualUnitaries])
        throw(ArgumentError("Invalid unitaries_type. Must be :Cliffords or :DualUnitaries"))
    end
    
    # Pre-compute bond lists (create once, use many times)
    # === ODD BONDS: [1,2], [3,4], [5,6], ... ===
    odd_bonds = []
    for i in 1:div(L, 2)
        push!(odd_bonds, (2*i-1, 2*i))
    end
    # Add periodic boundary condition for odd bonds if L is odd
    if is_pbc && (L % 2 == 1)
        push!(odd_bonds, (L, 1))
    end
    
    # === EVEN BONDS: [2,3], [4,5], [6,7], ... ===
    even_bonds = []
    for i in 1:div(L-1, 2)
        push!(even_bonds, (2*i, 2*i+1))
    end
    # Add periodic boundary condition for even bonds if L is even
    if is_pbc && (L % 2 == 0)
        push!(even_bonds, (L, 1))
    end
    
    # Track measurement outcomes: [time_step][qubit] -> outcome (or nothing if not measured)
    measurement_record = []
    gates_record = []
    
    # Main circuit loop
    for t in 1:T
        # === EVEN BONDS LAYER ===
        if extract_gates
            # If extracting gates, we need to store the gate indices
            gates_realization_even = Vector{Union{Nothing, Int}}(nothing, length(even_bonds))
        end
        # Apply random Clifford gates on even bonds
        for (bond_ind,(q1, q2)) in enumerate(even_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            if extract_gates; gates_realization_even[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_even); end
        
        # === FIRST MEASUREMENT LAYER ===
        measurement_outcomes_1 = Vector{Union{Nothing, Int}}(nothing, L)
        for qubit in 1:L
            if qubit == 1 && (!measure_first_qubit)
                continue  # Skip measurement of first qubit if not required
            elseif rand() < p_ar[qubit]  # Measure with probability p_x
                state, outcome, is_det = measure_z!(qubit, state)
                measurement_outcomes_1[qubit] = outcome
            end
        end
        push!(measurement_record, measurement_outcomes_1)

        # === ODD BONDS LAYER ===
        if extract_gates
            # If extracting gates, we need to store the gate indices
            gates_realization_odd = Vector{Union{Nothing, Int}}(nothing, length(odd_bonds))
        end
        # Apply random Clifford gates on odd bonds
        for (bond_ind,(q1, q2)) in enumerate(odd_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            if extract_gates; gates_realization_odd[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_odd); end
        
        # === SECOND MEASUREMENT LAYER ===
        measurement_outcomes_2 = Vector{Union{Nothing, Int}}(nothing, L)
        if t<T
            for qubit in 1:L
                if qubit == 1 && (!measure_first_qubit)
                    continue  # Skip measurement of first qubit if not required
                elseif rand() < p_ar[qubit]  # Measure with probability p_x
                    state, outcome, is_det = measure_z!(qubit, state)
                    measurement_outcomes_2[qubit] = outcome
                end
            end
        end
        push!(measurement_record, measurement_outcomes_2)
    end
    
    if extract_gates
        return state, measurement_record, gates_record
    else
        return state, measurement_record
    end
end


"""
Measures the entropy each time step.
"""
function brickwork_circuit_dynamics(L::Int, T::Int, p_ar::Vector{Float64}, state, is_pbc::Bool; extract_gates::Bool=false, measure_first_qubit::Bool=true, unitaries_type::Symbol=:Cliffords)
    # Validate inputs
    if !(length(p_ar)==L)
        throw(ArgumentError("probability vector p_ar must have length L"))
    end
    if !(unitaries_type in [:Cliffords, :DualUnitaries])
        throw(ArgumentError("Invalid unitaries_type. Must be :Cliffords or :DualUnitaries"))
    end
    
    # Pre-compute bond lists (create once, use many times)
    # === ODD BONDS: [1,2], [3,4], [5,6], ... ===
    odd_bonds = []
    for i in 1:div(L, 2)
        push!(odd_bonds, (2*i-1, 2*i))
    end
    # Add periodic boundary condition for odd bonds if L is odd
    if is_pbc && (L % 2 == 1)
        push!(odd_bonds, (L, 1))
    end
    
    # === EVEN BONDS: [2,3], [4,5], [6,7], ... ===
    even_bonds = []
    for i in 1:div(L-1, 2)
        push!(even_bonds, (2*i, 2*i+1))
    end
    # Add periodic boundary condition for even bonds if L is even
    if is_pbc && (L % 2 == 0)
        push!(even_bonds, (L, 1))
    end
    
    # Track measurement outcomes: [time_step][qubit] -> outcome (or nothing if not measured)
    measurement_record = []
    gates_record = []

    # Initialize system partitions for entropy calculations
    halfL = Vector{Int}(1:div(L, 2))  # First half of qubits
    A8 = Vector{Int}(1:div(L, 8))  # First antipodal eighth
    B8 = Vector{Int}(div(L, 2)+1:div(5*L, 8))  # Second antipodal eighth
    # Initialize observables record
    entropy_record = Vector{Int}(undef, 2*T)
    N_meas_record = Vector{Int}(undef, 2*T)
    N_det_record = Vector{Int}(undef, 2*T)
    I2_record = Vector{Int}(undef, 2*T)
    
    # Main circuit loop
    for t in 1:T
        # === EVEN BONDS LAYER ===
        if extract_gates
            # If extracting gates, we need to store the gate indices
            gates_realization_even = Vector{Union{Nothing, Int}}(nothing, length(even_bonds))
        end
        # Apply random Clifford gates on even bonds
        for (bond_ind,(q1, q2)) in enumerate(even_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            if extract_gates; gates_realization_even[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_even); end
        
        # === FIRST MEASUREMENT LAYER ===
        N_meas = 0 # Number of measurements in this layer
        N_det = 0  # Number of determinstic measurements in this layer 
        measurement_outcomes_1 = Vector{Union{Nothing, Int}}(nothing, L)
        for qubit in 1:L
            if qubit == 1 && (!measure_first_qubit)
                continue  # Skip measurement of first qubit if not required
            elseif rand() < p_ar[qubit]  # Measure with probability p_x
                state, outcome, is_det = measure_z!(qubit, state)
                measurement_outcomes_1[qubit] = outcome
                N_meas += 1
                if is_det; N_det += 1; end
            end
        end
        push!(measurement_record, measurement_outcomes_1)

        # Calculate entropy after first measurement layer
        entropy_record[2*t-1] = state_EE(state, halfL)
        N_meas_record[2*t-1] = N_meas
        N_det_record[2*t-1] = N_det
        I2_record[2*t-1] = state_EE(state, A8)+
                           state_EE(state, B8)-
                           state_EE(state, vcat(A8, B8))

        # === ODD BONDS LAYER ===
        if extract_gates
            # If extracting gates, we need to store the gate indices
            gates_realization_odd = Vector{Union{Nothing, Int}}(nothing, length(odd_bonds))
        end
        # Apply random Clifford gates on odd bonds
        for (bond_ind,(q1, q2)) in enumerate(odd_bonds)
            if unitaries_type == :Cliffords
                # Use Clifford gates
                gate_index = rand(0:11519)  # Random Clifford gate index
                apply_random_two_qubit_clifford!(q1, q2, state, gate_index=gate_index)
            else # unitaries_type == :DualUnitaries
                # Use dual-unitary gates (not implemented here, placeholder)
                gate_index = rand(0:5759)  # Placeholder for dual-unitary gate index
                apply_random_two_qubit_dual_unitary!(q1, q2, state, gate_index=gate_index)
            end
            if extract_gates; gates_realization_odd[bond_ind]=gate_index; end
        end
        if extract_gates; push!(gates_record, gates_realization_odd); end
        
        # === SECOND MEASUREMENT LAYER ===
        N_meas = 0 # Number of measurements in this layer
        N_det = 0  # Number of determinstic measurements in this layer
        measurement_outcomes_2 = Vector{Union{Nothing, Int}}(nothing, L)
        if t<T
            for qubit in 1:L
                if qubit == 1 && (!measure_first_qubit)
                    continue  # Skip measurement of first qubit if not required
                elseif rand() < p_ar[qubit]  # Measure with probability p_x
                    state, outcome, is_det = measure_z!(qubit, state)
                    measurement_outcomes_2[qubit] = outcome
                    N_meas += 1
                    if is_det; N_det += 1; end
                end
            end
        end
        push!(measurement_record, measurement_outcomes_2)

        # Calculate entropy after second measurement layer
        entropy_record[2*t] = state_EE(state, halfL)
        N_meas_record[2*t] = N_meas
        N_det_record[2*t] = N_det
        I2_record[2*t] = state_EE(state, A8)+
                         state_EE(state, B8)-
                         state_EE(state, vcat(A8, B8))
    end
    
    if extract_gates
        return state, entropy_record, N_meas_record, N_det_record, I2_record, measurement_record, gates_record
    else
        return state, entropy_record, N_meas_record, N_det_record, I2_record, measurement_record
    end
end



function calc_EE_observables(state::MixedDestabilizer)
    L = state.tab.nqubits
    halfL = Vector{Int}(1:div(L, 2))  # First half of qubits
    A_quarterL = Vector{Int}(1:div(L, 4))  # First quarter of qubits
    B_quarterL = Vector{Int}(div(L, 4)+1:div(L, 2))  # Second quarter of qubits
    C_quarterL = Vector{Int}(div(L, 2)+1:div(L, 2)+div(L, 4))  # Third quarter of qubits

    S = state_EE(state, halfL)
    I3 = state_EE(state, A_quarterL)+
         state_EE(state, B_quarterL)+
         state_EE(state, C_quarterL)-
         state_EE(state, vcat(A_quarterL, B_quarterL))-
         state_EE(state, vcat(B_quarterL, C_quarterL))-
         state_EE(state, vcat(A_quarterL, C_quarterL))+
         state_EE(state, vcat(A_quarterL, B_quarterL, C_quarterL))

    return S, I3
end

