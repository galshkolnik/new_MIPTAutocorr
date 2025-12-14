using new_MIPTAutocorr

L=80; 
T_therm = 4*L; T_meas=4*L
unitaries_type=:Cliffords; is_pbc=true
measure_I5=true

pmeas=0.16; initial_state_type=:product_rand
pmeas_ar = pmeas.*ones(L)

state_0 = new_MIPTAutocorr.create_initial_state(L, initial_state_type)

state_i, measurement_record = new_MIPTAutocorr.brickwork_circuit(L, T_therm, pmeas_ar, deepcopy(state_0), is_pbc; extract_gates=false, unitaries_type=unitaries_type)

state_f, entropy_record, N_meas_record, N_det_record, I2_record, I5_record, measurement_record = new_MIPTAutocorr.brickwork_circuit_dynamics(L, T_meas, pmeas_ar, deepcopy(state_i), is_pbc; extract_gates=false, unitaries_type=unitaries_type, measure_I5=measure_I5)

shannon_entropy = [sum(N_meas_record[1:t]) - sum(N_det_record[1:t]) for t in 1:length(N_meas_record)].*log(2)