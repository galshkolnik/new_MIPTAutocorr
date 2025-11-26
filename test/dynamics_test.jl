using new_MIPTAutocorr

L=32; 
T_therm = 2*L; T_meas=L
unitaries_type=:Cliffords; is_pbc=true

p=0.16; initial_state_type=:product_rand
pmeas_ar = p.*ones(L)

state_0 = new_MIPTAutocorr.create_initial_state(L, initial_state_type)

state_i, measurement_record = new_MIPTAutocorr.brickwork_circuit(L, T_therm, pmeas_ar, deepcopy(state_0), is_pbc; extract_gates=false, unitaries_type=unitaries_type)

state_f, entropy_record, N_meas_record, N_det_record, I2_record, measurement_record = new_MIPTAutocorr.brickwork_circuit_dynamics(L, T_meas, pmeas_ar, deepcopy(state_i), is_pbc; extract_gates=false, unitaries_type=unitaries_type)

shannon_entropy = [sum(N_meas_record[1:t]) - sum(N_det_record[1:t]) for t in 1:length(N_meas_record)].*log(2)

