#pragma once
#include <global.h>

ns_easyquantum

constexpr int pauliI = 0;
constexpr int pauliX = 1;
constexpr int pauliY = 2;
constexpr int pauliZ = 3;

constexpr int pauli2int(char p) {
	switch (p) {
	case 'I':
	case 'i':
		return pauliI;
	case 'X':
	case 'x':
		return pauliX;
	case 'Y':
	case 'y':
		return pauliY;
	case 'Z':
	case 'z':
		return pauliZ;
	default:
		assert(false, "Bad Pauli Operator.");
	}
	return -1;
}

template<typename qtraits>
struct Pauli {
	using qid = typename qtraits::qidx_t;
	using fp_t = typename qtraits::value_t; 
	using uint_t = typename qtraits::idx_t;
	using vector_t = DenseVector<fp_t, uint_t>;
	using matrix_t = DenseMatrix<fp_t, uint_t>;

	std::map<qid, char> paulis;

	Pauli() {}

	Pauli(std::map<qid, char> paulis_) : paulis(paulis_) { }
	Pauli(std::pair<qid, char> pauli) {
		paulis.insert(pauli);
	}

	Pauli(std::vector<std::pair<char, qid>> ps) {
		for (auto &p : ps) {
			auto iter = paulis.find(p.second);
			if (iter != paulis.end())
				assert(false, "Repeat qubits.");

			paulis.insert({ p.second, p.first });
		}
	}

	std::string to_string() const {
		if (paulis.size() == 0) return "";
		
		std::stringstream ss;
		
		for (auto iter = paulis.begin(); iter != std::prev(paulis.end()); ++iter) {
			ss << iter->second << iter->first << " ";
		}
		auto iter = std::prev(paulis.end());
		ss << iter->second << iter->first;
		return ss.str();
	}

	bool diagonal_only() const {
		for (auto& pauli : paulis) {
			int p = pauli2int(pauli.second);
			if (p == pauliX) return false;
			if (p == pauliY) return false;
			if (p == pauliZ || p == pauliI) continue;
			assert(false, "Bad pauli cases.");
		}
		return true;
	}

	bool check_space(std::vector<qid> trace_qubits) const {
		for (auto& pauli : paulis) {
			qid id = pauli.first;
			if (!is_in(trace_qubits, id)) {
				return false;
			}
		}
		return true;
	}

	static qid find_pos(const std::vector<qid>& trace_qubits, qid q) {
		auto iter = std::find(trace_qubits.begin(), trace_qubits.end(), q);
		assert(iter != trace_qubits.end(), 
			"Bad calling sequence. Call this after check_space();");

		return qid(iter - trace_qubits.begin());
	}

	void append_diagonal(std::vector<qid> trace_qubits, fp_t coef, vector_t& vec) const {
		uint_t size = pow2<uint_t>(trace_qubits.size());
		assert(vec.size == size, "Bad calling sequence.");
		for (vector_t::size_type i = 0; i < size; ++i) {
			int parity = 0; // even case
			for (auto pauli : paulis) {
				qid pos = find_pos(trace_qubits, pauli.first);
				int p = pauli2int(pauli.second);
				assert(p == pauliI || p == pauliZ, "Bad calling sequence. Call this after diagonal_only()");
				if (p == pauliI) continue;
				// get parity of i at pos
				parity += (int(i >> pos));
			}
			parity &= 1;
			if (parity == 0) vec[i] += coef;
			if (parity == 1) vec[i] -= coef;
		}
	}

	void append_non_diagonal(std::vector<qid> trace_qubits, fp_t coef, vector_t& vec) const {
		uint_t size = pow2<uint_t>(trace_qubits.size());
		assert(vec.size == size, "Bad calling sequence.");
		for (vector_t::size_type i = 0; i < size; ++i) {
			int parity = 0; // even case
			for (auto pauli : paulis) {
				qid pos = find_pos(trace_qubits, pauli.first);
				int p = pauli2int(pauli.second);
				// assert(p == pauliI || p == pauliZ, "Bad calling sequence. Call this after diagonal_only()");
				if (p == pauliI) continue;
				// get parity of i at pos
				parity += (int(i >> pos));
			}
			parity &= 1;
			if (parity == 0) vec[i] += coef;
			if (parity == 1) vec[i] -= coef;
		}
	}

	/* Warning: Read this before use this function!! 
	   
	   When we get a measurement result and calculate the expecatation,
	   the pauliX/pauliY term will need to append a Hadamard/X(pi/2) to
	   the end of the circuit. This function assumes your measurement
	   result is from the circuit which the gate is appended to.
	*/
	int get_expectation_from_measurement(size_t meas) const {
		int res = 1;
		for (auto& pauli : paulis) {
			if (res == 0) return 0;
			int p = pauli2int(pauli.second);
			qid q = pauli.first;			
			bool value = (meas >> q) % 2;
			switch (p) {
			case pauliX:
			case pauliY:
			case pauliZ:
				if (value) res *= -1;
				break;
			case pauliI:
				break;
			default:
				assert(false);
			}
		}
		return res;
	}
};

template<typename qtraits = default_qtraits>
void circuit_append(Circuit<qtraits>& c, const Pauli<qtraits>& paulis) {
	for (const auto& pauli : paulis.paulis) {
		int p = pauli2int(pauli.second);
		switch (p) {
		case pauliX:
			c - H(pauli.first);
			break;
		case pauliY:
			c - RX(pauli.first, pi / 2);
			break;
		case pauliZ:
		case pauliI:
			break;
		default:
			assert(false);
		}
	}
}

template<typename qtraits>
struct QuantumLess {
	using qid = typename qtraits::qidx_t;
	
	/* I < X < Y < Z */
	constexpr bool operator()(const char p1, const char p2) const {
		return pauli2int(p1) < pauli2int(p2);
	}

	constexpr bool operator()(const std::map<qid, char>& p1, 
		const std::map<qid, char>& p2) const {

		size_t size1 = p1.size();
		size_t size2 = p2.size();

		auto iter1 = p1.begin();
		auto iter2 = p2.begin();
		for (size_t i = 0;; ++i) {
			if (iter1 == p1.end() && iter2 == p2.end()) return false;
			if (iter1 == p1.end()) return true;
			if (iter2 == p2.end()) return false;
			if (iter1->first < iter2->first) return true;
			if (iter1->first > iter2->first) return false;
			if (iter1->second < iter2->second) return true;
			if (iter1->second > iter2->second) return false;

			iter1++; iter2++;
		}
		return false;
	}

	constexpr bool operator()(const Pauli<qtraits>& p1, const Pauli<qtraits>& p2) const {
		return p1.paulis < p2.paulis;
	}

	constexpr bool operator()(const std::vector<std::pair<char, qid>>& p1,
		const std::vector<std::pair<char, qid>>& p2) const {
		Pauli<qtraits> pauli1(p1);
		Pauli<qtraits> pauli2(p2);
		return pauli1.paulis < pauli2.paulis;
	}
};

template<typename qtraits = default_qtraits>
struct Operator {
	using fp_t = typename qtraits::value_t;
	using complex_t = std::complex<fp_t>;
	using pauli_t = Pauli<qtraits>;
	using op_terms_t = std::map<Pauli<qtraits>, complex_t, QuantumLess<qtraits>>;
	using op_term_t = std::pair<Pauli<qtraits>, complex_t>;
	using ham_terms_t = std::map<Pauli<qtraits>, fp_t, QuantumLess<qtraits>>;
	using ham_term_t = std::pair<Pauli<qtraits>, fp_t>;
	op_terms_t terms;

	Operator() {}

	Operator(const pauli_t& pauli) {
		terms.insert({ pauli, 1 });
	}

	Operator(const op_term_t& term) {
		terms.insert(term);
	}

	Operator(const op_terms_t& terms_) {
		terms = terms_;
	}

	Operator(const Operator<qtraits>& op) {
		terms = op.terms;
	}

	complex_t& operator[](const pauli_t& pauli) {
		return terms[pauli];
	}

	void clear() {
		auto iter = std::remove_if(terms.begin(), terms.end(),
			[](const op_term_t& op) {
				if (round0(op.second, quantum_default_threshold))
					return true;
				return false;
			});

		terms.erase(iter, terms.end());
	}

	static ham_terms_t to_hamiltonian(const op_terms_t& terms) {
		ham_terms_t ret;
		for (auto& term : terms) {
			if (round0(term.second.imag, default_quantum_threshold)) {
				ret[term.first] = term.second.real;
			}
			else assert(false, "There exists an imaginary coefficient.");
		}
		return ret;
	}

	Operator<qtraits>& operator+=(const Operator<qtraits>& op) {
		Operator<qtraits> &newop = *this;
		for (const auto& term : op.terms) {
			auto iter = newop.terms.find(term.first);
			if (iter == newop.terms.end()) {
				newop.terms[term.first] = 0;
			}
			newop.terms[term.first] += term.second;
		}
		return newop;
	}

	Operator<qtraits> operator+(const Operator<qtraits>& op) {
		Operator<qtraits> newop(*this);
		newop += op;
		return newop;
	}

	Operator<qtraits>& operator-=(const Operator<qtraits>& op) {
		Operator<qtraits>& newop = *this;
		for (const auto& term : op.terms) {
			auto iter = newop.terms.find(term.first);
			if (iter == newop.terms.end()) {
				newop.terms[term.first] = 0;
			}
			newop.terms[term.first] -= term.second;
		}
		return newop;
	}

	Operator<qtraits> operator-(const Operator<qtraits>& op) {
		Operator<qtraits> newop(*this);
		newop -= op;
		return newop;
	}

	Operator<qtraits>& operator*=(complex_t c) {
		Operator<qtraits>& newop = *this;
		for (auto& term : newop.terms) {
			term.second *= c;
		}
		return newop;
	}

	Operator<qtraits> operator*(complex_t c) {
		Operator<qtraits> newop(*this);
		newop *= c;
		return newop;
	}

	std::string to_string() const {
		std::function<std::string(const std::pair<Pauli<qtraits>, complex_t>&)> func =
			[](const std::pair<Pauli<qtraits>, complex_t>& term)
		{return complex_to_string(term.second) + " " + term.first.to_string(); };
		return vec2str(terms, func, std::string(" + "));
	}
};

template<typename qtraits>
struct Hamiltonian {
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;
	using uint_t = typename qtraits::idx_t;
	using vector_t = DenseVector<fp_t, uint_t>;
	using matrix_t = DenseMatrix<fp_t, uint_t>;

	vector_t full_diagonal;
	// full_diagonal_valid > 0  : full_diagonal
	// full_diagonal_valid == 0 : never_check
	// full_diagonal_valid < 0  : non diagonal
	int full_diagonal_valid = 0;

	std::map<Pauli<qtraits>, fp_t, QuantumLess<qtraits>> terms;

	Hamiltonian(Operator<qtraits> &op) {
		terms = Operator<qtraits>::to_hamiltonian(op.terms);
	}

	Hamiltonian(std::vector<std::pair<char, qid>> ps) {
		terms.insert({ ps, 1. });
	}

	Hamiltonian(fp_t coef, std::vector<std::pair<char, qid>> ps) {
		terms.insert({ ps, coef });
	}

	Hamiltonian(std::map<std::vector<std::pair<char, qid>>, fp_t, QuantumLess<qtraits>> ps) {
		for (auto& p : ps) {
			bool is_success = terms.insert({ p.first, p.second }).second;
			assert(is_success, "Duplicate Key");
		}
	}

	Hamiltonian(const Hamiltonian& h) {
		if (h.full_diagonal_valid > 0) {
			full_diagonal = h.full_diagonal;
			full_diagonal_valid = h.full_diagonal_valid;
		}
		terms = h.terms;
	}

	bool get_full_diagonal(std::vector<qid> trace_qubits) {
		if (full_diagonal_valid)
			return full_diagonal_valid > 0;

		uint_t size = pow2<uint_t>(trace_qubits.size());

		bool first_create = true;
		full_diagonal_valid = 1;
		for (const auto& term : terms) {
			if (!term.first.diagonal_only()) {				
				full_diagonal_valid = -1;
				break;
			}
		}
		full_diagonal.initialize(size);

		if (full_diagonal_valid > 0) {
			for (const auto& term : terms) {
				term.first.append_diagonal(trace_qubits, term.second, full_diagonal);
			}
		}

		return full_diagonal_valid > 0;
	}

	/*bool get_full_matrix() {
		return full_matrix_valid;
	}*/

	std::string to_string() const {
		std::function<std::string(const std::pair<Pauli<qtraits>, fp_t>&)> func =
			[](const std::pair<Pauli<qtraits>, fp_t>& term) 
		{return std::to_string(term.second) + " " + term.first.to_string(); };
		return vec2str(terms, func, std::string(" + "));
	}

	/* Only available when all Z*/
	fp_t get_expectation_from_measurement(size_t meas_res) const {
		// size_t nterms = terms.size();
		assert(full_diagonal_valid > 0, "Only diagonal is valid.");
		fp_t value = 0;
		for (auto& term : terms) {
			fp_t energy = term.second * term.first.get_expectation_from_measurement(meas_res);
			value += energy;
			std::cout << "meas: " << meas_res << " value: " << energy << std::endl;
		}
		return value;
	}

	fp_t get_expectation_from_quantum_state(fp_t* real, fp_t* imag, uint_t size, const std::vector<qid> &trace_qubits) {
		assert(full_diagonal_valid < 0, "Only non-diagonal is valid.");
		fp_t val = 0;
		for (const auto& term : terms) {
			for (const auto& pauli : term.first.paulis) {
				int p = pauli2int(pauli.second);
				switch (p) {
				case pauliX:
					state_manipulator<qtraits>::h(real, imag, size, pauli.first);
					break;
				case pauliY:
					state_manipulator<qtraits>::rx(real, imag, size, pi / 2, pauli.first);
					break;
				case pauliZ:
				case pauliI:
					break;
				default:
					assert(false);
				}
			}			
			term.first.append_non_diagonal(trace_qubits, term.second, full_diagonal);
			val += result_analyzer<qtraits>::get_expectation_from_amplitude(
				real, imag, size, full_diagonal.data);

			for (const auto& pauli : term.first.paulis) {
				int p = pauli2int(pauli.second);
				switch (p) {
				case pauliX:
					state_manipulator<qtraits>::h(real, imag, size, pauli.first);
					break;
				case pauliY:
					state_manipulator<qtraits>::rx(real, imag, size, -pi / 2, pauli.first);
					break;
				case pauliZ:
				case pauliI:
					break;
				default:
					assert(false);
				}
			}

			full_diagonal.clear();
		}
		return val;
	}

	fp_t get_expectation_from_circuit_montecarlo_singlethread(
		const Circuit<qtraits>& circuit, size_t shots, RandomEngine* rng) {

		assert(full_diagonal_valid < 0, "Use get_expectation_from_measurement instead.");

		fp_t exp = 0;
		for (auto& term : terms) {
			// copy the circuit
			Circuit<qtraits> c(circuit);
			// append according to the Pauli term
			int val = 0;
			circuit_append(c, term.first);
			uint_t* res = simulator_v1<qtraits>::simulate_N(shots, c, rng);

			for (size_t i = 0; i < shots; ++i) {
				val += term.first.get_expectation_from_measurement(res[i]);
			}
			exp += term.second * val;
		}
		exp /= shots;
		return exp;
	}

	fp_t get_expectation_from_circuit_montecarlo_multithread(
		const Circuit<qtraits>& circuit, size_t shots, RandomEngine* rng) {

		assert(full_diagonal_valid < 0, "Use get_expectation_from_measurement instead.");

		fp_t exp = 0;
		for (auto& term : terms) {
			// copy the circuit
			Circuit<qtraits> c(circuit);
			// append according to the Pauli term
			int val = 0;
			circuit_append(c, term.first);
			uint_t* res = simulator_v1<qtraits>::simulate_N_threads(shots, c, rng);

			for (size_t i = 0; i < shots; ++i) {
				val += term.first.get_expectation_from_measurement(res[i]);
			}
			exp += term.second * val;
		}
		exp /= shots;
		return exp;
	}
};



ns_end
