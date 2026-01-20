using new_MIPTAutocorr

L=32; 
T_therm = 4*L; T_meas=4*L
unitaries_type=:Cliffords; is_pbc=true
measure_I5=false; measure_I3=true

pmeas=0.14; initial_state_type=:product_rand
pmeas_ar = pmeas.*ones(L)

state_0 = new_MIPTAutocorr.create_initial_state(L, initial_state_type)

state_i, measurement_record = new_MIPTAutocorr.brickwork_circuit(L, T_therm, pmeas_ar, deepcopy(state_0), is_pbc; unitaries_type=unitaries_type)

res = new_MIPTAutocorr.brickwork_circuit_dynamics(L, T_meas, pmeas_ar, deepcopy(state_i), is_pbc; unitaries_type=unitaries_type, measure_I5=measure_I5, measure_I3=measure_I3)

state_f = res.state
entropy_record = res.entropy_record
N_meas_record = res.N_meas_record
N_det_record = res.N_det_record
I3_record = res.I3_record