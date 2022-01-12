#pragma once

#include <global.h>
#include <quantum.h>
#include "Hamiltonian.h"

ns_easyquantum

template<typename qtraits_t = default_qtraits>
struct ObjectiveFunction {
	using qtraits = qtraits_t;
	using fp_t = typename qtraits::value_t;
	using uint_t = typename qtraits::idx_t;
	using qid = typename qtraits::qidx_t;
	std::vector<Variable<qtraits>> variables;
	Circuit<qtraits_t> circuit;
	RealCircuit<qtraits_t> real_circuit;
	Hamiltonian<qtraits_t> hamiltonian;

	RandomEngine** rngs = nullptr;
	uint_t* res = nullptr;
	fp_t* real = nullptr;
	fp_t* imag = nullptr;
	size_t shots = 0;

	qid qn = 0;
	uint_t state_size = 0;

	ObjectiveFunction(std::vector<Variable<qtraits_t>> vars,
		Circuit<qtraits_t> circ, Hamiltonian<qtraits_t> hami)
		: variables(vars), circuit(circ), 
		hamiltonian(hami) {}

	Circuit<qtraits>* get_circuit_ptr() {
		return std::addressof(circuit);
	}

	void ready() {
		real_circuit.c = circuit.valuecopy();
		real_circuit.ready();
	}

	void refresh_all_variables() {
		for (auto& var : variables) var.refresh(get_circuit_ptr());
	}

	void create_space(size_t shots, qid n_qubit) {
		if (this->shots != 0) {
			assert(shots == this->shots, "Cannot perform 2 calculation with different shots. "
				"Use another ObjectiveFunction instead.");
		}
		this->shots = shots; 
		if (rngs == nullptr) {
			rngs = make_engines(shots);
		}

		if (real == nullptr) {
			// When first calling this function, we create the space for real, imag and res.
			assert(imag == nullptr, "Why real == nullptr and imag != nullptr? Weird case.");
			assert(res == nullptr, "Why real == nullptr and res != nullptr? Weird case.");
			res = new uint_t[shots];

			assert(qn == 0 && state_size == 0, "Created state before? Weird case.");
			qn = n_qubit;
			state_size = pow2<uint_t, qid>(qn);

			real = new fp_t[shots * state_size];
			imag = new fp_t[shots * state_size];
		}
		else {
			// If it is not the first call, we will not create.
			assert(imag != nullptr, "Why real != nullptr and imag == nullptr? Weird case.");
			assert(res != nullptr, "Why real == nullptr and res == nullptr? Weird case.");

			assert(qn == n_qubit && state_size == pow2<uint_t, qid>(n_qubit),
				"Cannot perform 2 calculation with different qubits. "
				"Use another ObjectiveFunction instead.");
		}
	}

	void simulate_shots(const Circuit<qtraits>& c, const RealCircuit<qtraits>& rc,
		bool thread, bool noisy) {
		if (noisy) {
			if (thread)
				simulator_v1<qtraits>::simulate_N_threads_noisy(shots, rc, rngs, real, imag, res);
			else
				simulator_v1<qtraits>::simulate_N_noisy(shots, rc, rngs, real, imag, res);
		}
		else {
			if (thread)
				simulator_v1<qtraits>::simulate_N_threads(shots, c, rngs, real, imag, res);
			else
				simulator_v1<qtraits>::simulate_N(shots, c, rngs, real, imag, res);
		}
	}

	fp_t _get_expectation_montecarlo(size_t shots, bool thread, bool noise) {
		bool diagonal_mode = hamiltonian.get_full_diagonal(
			result_analyzer<qtraits>::generate_qubit_idxes(circuit.max_qubit));

		create_space(shots, circuit.max_qubit);
		
		if (diagonal_mode) {
			if (noise)
				real_circuit.ready();
			simulate_shots(circuit, real_circuit, thread, noise);
			fp_t sum = 0;
			for (size_t i = 0; i < shots; ++i) {
				sum += hamiltonian.full_diagonal[res[i]];
			}
			sum /= shots;
			return sum;
		}
		else {
			fp_t exp = 0;
			for (const auto& term : hamiltonian.terms) {
				Circuit<qtraits> c = circuit.valuecopy();
				circuit_append(c, term.first);
				RealCircuit<qtraits> rc(c);

				if (noise) {
					rc.ready();
				}
				
				int val = 0;
				simulate_shots(c, rc, thread, noise);
				for (size_t i = 0; i < shots; ++i) {
					val += term.first.get_expectation_from_measurement(res[i]);
				}
				exp += term.second * val;
			}
			exp /= shots;
			return exp;
		}
	}

	fp_t _get_expectation_prob_noisefree() {
		std::vector<qid> qubits = 
			result_analyzer<qtraits>::generate_qubit_idxes(circuit.max_qubit);

		bool diagonal_mode = hamiltonian.get_full_diagonal(qubits);

		fp_t* real = nullptr, *imag = nullptr;

		uint_t size = state_manipulator<qtraits>::create_state(
			&real, &imag, circuit.max_qubit);

		simulator_v1<qtraits>::run_circuit(circuit, real, imag, size);

		fp_t res = 0;

		if (diagonal_mode) {
			res = result_analyzer<qtraits>::get_expectation_from_amplitude(
				real, imag, size, hamiltonian.full_diagonal.data);
		}
		else {
			res = hamiltonian.get_expectation_from_quantum_state(
				real, imag, size, qubits);
		}
		free(real);
		free(imag);
		return res;
	}

	fp_t get_expectation_montecarlo(size_t shots, bool thread, bool noise) {
		refresh_all_variables();
		return _get_expectation_montecarlo(shots, thread, noise);
	}

	fp_t get_derivative_montecarlo(size_t n_var, size_t shots, bool thread, bool noise) {
		fp_t derivative = 0;
		// for i-th variable
		refresh_all_variables();

		for (size_t i = 0; i < variables[n_var].size(get_circuit_ptr()); ++i) {
			variables[n_var].temporarily_set(get_circuit_ptr(), i, pi / 2);
			derivative += (_get_expectation_montecarlo(shots, thread, noise) / 2);
			variables[n_var].temporarily_set(get_circuit_ptr(), i, -pi / 2);
			derivative -= (_get_expectation_montecarlo(shots, thread, noise) / 2);
			variables[n_var].refresh(get_circuit_ptr());
		}
		return derivative;
	}

	std::vector<fp_t> get_derivatives_montecarlo(size_t shots, bool thread, bool noise) {
		std::vector<fp_t> derivatives(variables.size());
		for (size_t i = 0; i < variables.size(); ++i) {
			derivatives[i] = get_derivative_montecarlo(i, shots, thread, noise);
		}
		return derivatives;
	}

	fp_t get_expectation_prob_noisefree() {
		refresh_all_variables();
		return _get_expectation_prob_noisefree();
	}

	fp_t get_derivative_prob_noisefree(size_t n_var) {
		fp_t derivative = 0;
		// for i-th variable
		refresh_all_variables();

		for (size_t i = 0; i < variables[n_var].size(get_circuit_ptr()); ++i) {
			variables[n_var].temporarily_set(get_circuit_ptr(), i, pi / 2);
			derivative += (_get_expectation_prob_noisefree() / 2);
			variables[n_var].temporarily_set(get_circuit_ptr(), i, - pi / 2);
			derivative -= (_get_expectation_prob_noisefree() / 2);
			variables[n_var].refresh(get_circuit_ptr());
		}
		return derivative;
	}

	std::vector<fp_t> get_derivatives_prob_noisefree() {
		std::vector<fp_t> derivatives(variables.size());
		for (size_t i = 0; i < variables.size(); ++i) {
			derivatives[i] = get_derivative_prob_noisefree(i);
		}
		return derivatives;
	}

	~ObjectiveFunction() {
		if (rngs != nullptr)
			free_engines(rngs, shots);
		if (res != nullptr)
			delete(res);
		if (real != nullptr)
			delete(real);
		if (imag != nullptr)
			delete(imag);
	}

};

ns_end