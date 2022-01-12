#pragma once

#include <global.h>
#include <randomutil.h>
#include <mathutil.h>
#include <testutil.h>
#include <iomanip>
#include "state_manipulator.h"

ns_easyquantum

#define TO_STRING_CASE(gatetype) case GateType::gatetype:
#define TO_STRING_DAG if (dag) {ss<<"dag ";}

/* gate type definition */
enum class GateType {
	null,
	CNOT,
	RX,
	RY,
	RZ,
	H,
	DIAGONAL,
	I,
};

constexpr double quantum_default_threshold = 1e-8;

// forward decl
template<typename T> class Circuit;
template<typename T> class Gate;

template<typename qtraits_t = default_qtraits>
struct VariableImpl {
	using qtraits = qtraits_t;
	using circuit_t = Circuit<qtraits>;
	using gate_t = Gate<qtraits>;
	using fp_t = typename qtraits::value_t;
	using gaterefs_t = std::vector<Gate<qtraits_t>*>;
	using tape_t = std::map<Circuit<qtraits_t>*, gaterefs_t>;
	using tape_iterator_t = typename tape_t::iterator;

	tape_t tape;
	fp_t value;

	VariableImpl() { value = 0; }
	VariableImpl(double setval) { value = setval; }
	VariableImpl(const VariableImpl<qtraits>& v) = delete;

	tape_iterator_t find_circuit(circuit_t* c) {
		auto citer = tape.find(c);
		assert(citer != tape.end(), "Trying to find a circuit not recorded.");
		return citer;
	}

	size_t size(circuit_t* c) {	
		find_circuit(c);
		return tape[c].size();
	}

	void add_ref(circuit_t* c, gate_t* g) {
		tape[c].push_back(g);
	}

	void dec_ref(circuit_t* c) {
		auto citer = tape.find(c);
		if (citer != tape.end())
			tape.erase(citer);
	}

	void refresh(circuit_t* c) {
		auto citer = find_circuit(c);
		for (auto gate : tape[c]) {
			gate->argument = value;
		}
	}

	void temporarily_set(circuit_t *c, size_t n_th, fp_t offset) {		
		for (size_t i = 0; i < size(c); ++i) {
			if (i != n_th)
				tape[c][i]->argument = value;
			else
				tape[c][i]->argument = value + offset;
		}
	}

	void set_value(fp_t new_value) {
		value = new_value;
	}
};

template<typename qtraits_t = default_qtraits>
struct Variable {
	using fp_t = typename qtraits_t::value_t;
	using circuit_t = typename VariableImpl<qtraits_t>::circuit_t;
	std::shared_ptr<VariableImpl<qtraits_t>> impl;

	Variable(const Variable<qtraits_t>& v) {
		impl = v.impl;
	}
	Variable(fp_t value = 0) {
		impl = std::make_shared<VariableImpl<qtraits_t>>(value);
	}
	VariableImpl<qtraits_t> * operator->() {
		return impl.get();
	}
	void set_value(fp_t value) {
		impl->set_value(value);
	}
	void set_offset(fp_t offset) {
		impl->set_value(impl->value + offset);
	}
	size_t size(circuit_t* c) {
		return impl->size(c);
	}
	void temporarily_set(circuit_t* c, size_t n_th, fp_t offset) {
		impl->temporarily_set(c, n_th, offset);
	}
	void refresh(circuit_t* c) {
		impl->refresh(c);
	}
};

template<typename qtraits = default_qtraits>
class Gate {
public:
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;

	std::vector<qid> qubits;
	fp_t argument;
	Variable<qtraits> var;
	bool isvar = false;
	GateType type = GateType::null;
	std::vector<qid> controllers;
	bool dag = false;

	Gate() {}

	Gate(GateType type_, std::vector<qid> qs, fp_t arg = 0)
		: type(type_),
		qubits(qs.begin(), qs.end()),
		argument(arg) {}

	Gate(GateType type_, std::vector<qid> qs, Variable<qtraits> var_)
		: type(type_),
		qubits(qs.begin(), qs.end()),
		var(var_),
		argument(var_->value),
		isvar(true) {}

	Gate(const Gate<qtraits>& g) {
		qubits.assign(g.qubits.begin(), g.qubits.end());
		argument = g.argument;
		type = g.type;
		controllers.assign(g.controllers.begin(), g.controllers.end());
		dag = g.dag;
		var = g.var;
		isvar = g.isvar;
	}

	Gate<qtraits> valuecopy() const {
		Gate<qtraits> g;
		g.qubits = qubits;
		g.argument = argument;
		g.type = type;
		g.controllers = controllers;
		g.dag = dag;
		return g;
	}

	Gate<qtraits>& record(Variable<qtraits>& v) {
		var = v;
		isvar = true;
		return *this;
	}

	void dagger() {
		dag = !dag;
	}

	Gate<qtraits> get_dagger() const {
		Gate<qtraits> g(*this);
		g.dag = !g.dag;
		return g;
	}

	void control(std::vector<qid> new_controllers) {
		controllers.insert(g.controllers.end(), 
			new_controllers.begin(), new_controllers.end());
	}

	Gate<qtraits> get_control(std::vector<qid> new_controllers) const {
		Gate<qtraits> g(*this);
		g.controllers.insert(g.controllers.end(), new_controllers.begin(), new_controllers.end());
		return g;
	}

	static std::string qasm_style_qubit(qid i) {
		std::stringstream ss;
		ss << "q[" << i << "]";
		return ss.str();
	}
	static std::string controller_to_string(std::vector<qid> controllers) {
		if (controllers.size() == 0) return std::string();
		std::stringstream ss;
		ss << "ctrl:{" << vec2str(controllers) << "}";
		return ss.str();
	}
	static std::string params_to_string(std::vector<fp_t> args) {
		if (args.size() == 0) return std::string();
		std::stringstream ss;
		ss << std::showpoint;
		ss << "(" << vec2str(args, ", ", SHOWPOINT) << ")";
		return ss.str();
	}

	std::string to_string() const {
		std::stringstream ss;
		switch (type) {
			TO_STRING_CASE(CNOT) {
				ss << "CNOT ";
				ss << qubits[0] << ", " << qubits[1] << " ";
				ss << controller_to_string(controllers);
				return ss.str();
			}
			TO_STRING_CASE(RX) {
				ss << "RX ";
				TO_STRING_DAG;
				ss << qubits[0] << " ";
				ss << params_to_string({ argument });
				ss << " ";
				ss << controller_to_string(controllers);
				return ss.str();
			}
			TO_STRING_CASE(RY) {
				ss << "RY ";
				TO_STRING_DAG;
				ss << qubits[0] << " ";
				ss << params_to_string({ argument });
				ss << " ";
				ss << controller_to_string(controllers);
				return ss.str();
			}
			TO_STRING_CASE(RZ) {
				ss << "RZ ";
				TO_STRING_DAG;
				ss << qubits[0] << " ";
				ss << params_to_string({ argument });
				ss << " ";
				ss << controller_to_string(controllers);
				return ss.str();
			}
			TO_STRING_CASE(H) {
				ss << "H ";
				ss << qubits[0] << " ";
				ss << controller_to_string(controllers);
				return ss.str();
			}
			TO_STRING_CASE(DIAGONAL) {
				ss << "Diagonal ";
				ss << 0 << ";";
				return ss.str();
			}
			TO_STRING_CASE(I) {
				ss << "I ";
				ss << qubits[0] << " ";
				ss << params_to_string({ argument });
				return ss.str();
			}
		default: {
			assert(false, "Bad Type.");
		}
		}
	}
	std::string to_qasm() const {
		assert(controllers.size() == 0);
		assert(dag == false);

		std::stringstream ss;
		switch (type) {
			TO_STRING_CASE(CNOT) {
				ss << "cx ";
				ss << qasm_style_qubit(qubits[0]) << ", " << qasm_style_qubit(qubits[1]);
				return ss.str();
			}
			TO_STRING_CASE(RX) {
				ss << "rx ";
				ss << params_to_string({ argument }) << " ";
				ss << qasm_style_qubit(qubits[0]);
				return ss.str();
			}
			TO_STRING_CASE(RY) {
				ss << "ry ";
				ss << params_to_string({ argument }) << " ";
				ss << qasm_style_qubit(qubits[0]);
				return ss.str();
			}
			TO_STRING_CASE(RZ) {
				ss << "rz ";
				ss << params_to_string({ argument }) << " ";
				ss << qasm_style_qubit(qubits[0]);
				return ss.str();
			}
			TO_STRING_CASE(H) {
				ss << "h ";
				ss << qasm_style_qubit(qubits[0]);
				return ss.str();
			}
			TO_STRING_CASE(DIAGONAL) {
				assert(false, "Bad Type.");
			}
			TO_STRING_CASE(I) {
				ss << "Id ";
				ss << qubits[0] << " ";
				return ss.str();
			}
			default: assert(false, "Bad Type."); 
		}
		return "";
	}

	~Gate() { }
};

template<typename qtraits = default_qtraits>
Gate<qtraits> CNOT(typename qtraits::qidx_t controller, typename qtraits::qidx_t target) {
	return Gate(GateType::CNOT, { controller, target });
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RX(typename qtraits::qidx_t q, typename qtraits::value_t arg = 0) {
	return Gate<qtraits>(GateType::RX, { q }, arg);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RY(typename qtraits::qidx_t q, typename qtraits::value_t arg) {
	return Gate(GateType::RY, { q }, arg);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RZ(typename qtraits::qidx_t q, typename qtraits::value_t arg) {
	return Gate(GateType::RZ, { q }, arg);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> H(typename qtraits::qidx_t q) {
	return Gate(GateType::H, { q });
}

template<typename qtraits = default_qtraits>
Gate<qtraits> I(typename qtraits::qidx_t q, size_t time) {
	using fp_t = typename qtraits::value_t;
	return Gate(GateType::I, { q }, (fp_t)time);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RX(typename qtraits::qidx_t q, Variable<qtraits> arg) {
	return Gate<qtraits>(GateType::RX, { q }, arg);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RY(typename qtraits::qidx_t q, Variable<qtraits> arg) {
	return Gate<qtraits>(GateType::RY, { q }, arg);
}

template<typename qtraits = default_qtraits>
Gate<qtraits> RZ(typename qtraits::qidx_t q, Variable<qtraits> arg) {
	return Gate<qtraits>(GateType::RZ, { q }, arg);
}

class flatten {};

template<typename qtraits = default_qtraits>
class Circuit {
public:
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;
	
	qid max_qubit = 0;
	std::list<Gate<qtraits>> gates;

	Circuit() {}

	void _circuit_append(const Circuit<qtraits>& c) {
		if (c.max_qubit > max_qubit)
			max_qubit = c.max_qubit;
		for (auto gate : c.gates) {
			gates.push_back(gate);
			if (gate.isvar) {
				gate.var->add_ref(this, std::addressof(gates.back()));
			}
		}
	}

	Circuit(const Circuit<qtraits>& c) {
		_circuit_append(c);
	}

	Circuit valuecopy() const {
		Circuit c;
		c.max_qubit = max_qubit;
		for (const auto& g : gates) {
			c.gates.push_back(g.valuecopy());
		}
		return c;
	}

	Circuit& operator-(Gate<qtraits> g) {
		if (g.type != GateType::DIAGONAL) {
			for (auto qubit : g.qubits) {
				if (max_qubit - 1 < qubit)
					max_qubit = qubit + 1;
			}
		}
		gates.push_back(g);
		if (g.isvar) {
			g.var->add_ref(this, std::addressof(gates.back()));
		}
		return *this;
	}

	Circuit& operator-(flatten flatten) {
		for (auto &gate : gates) {
			assert(!gate.isvar, "Not support variable gate flatten.");
			if (gate.dag) {
				gate.dag = false;
				switch (gate.type) {
				case GateType::RX:
				case GateType::RY:
				case GateType::RZ:
					gate.argument *= -1;
					break;
				default:
					break;
				}
			}
		}
		return *this;
	}

	Circuit& operator-(const Circuit<qtraits> &c) {
		_circuit_append(c);
		return *this;
	}

	Circuit& dagger() {
		gates.reverse();
		for (auto& gate : gates) {
			gate.dagger();
		}
		return newc;
	}

	Circuit control(std::vector<qid> controllers) {
		Circuit newc(*this);
		for (auto& gate : newc.gates) {
			gate.control(controllers);
		}
		return newc;
	}

	std::string to_string() const {
		std::stringstream ss_stat;
		std::stringstream ss_header;
		std::stringstream ss_groupdefs;

		/* header*/
		ss_header << "qubit : " << max_qubit << ";" << std::endl;

		/* stat */
		for (auto gate : gates) {
			ss_stat << gate.to_string() << ";" << std::endl;
		}
		return ss_header.str() + "\n" + ss_stat.str() + "\n";
	}

	std::string to_qasm() const {
		std::stringstream ss;
		ss << "OPENQASM 2.0;" << std::endl
			<< "include \"qelib1.inc\";" << std::endl
			<< "qreg q[" << max_qubit << "];" << std::endl;

		for (const auto &gate : gates) {
			ss << gate.to_qasm() << ";" << std::endl;
		}
		return ss.str();
	}

	~Circuit() {
		for (auto& gate : gates) {
			if (gate.isvar) {
				gate.var->dec_ref(this);
			}
		}
	}
};

void get_damping_kraus_op(
	double* k0_real, double* k0_imag,
	double* k1_real, double* k1_imag,
	const int T1, const int T_gate);

void get_dephasing_kraus_op(
	double* k0_real, double* k0_imag,
	double* k1_real, double* k1_imag,
	const int T1, const int T2, const int T_gate);

template<typename qtraits = default_qtraits>
class RealCircuit {
public:
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;

	struct Kraus {
		fp_t *kraus0_real;
		fp_t *kraus0_imag;
		fp_t *kraus1_real;
		fp_t *kraus1_imag;
	};
	std::vector<Kraus> one_qubit_damping_kraus;
	std::vector<Kraus> one_qubit_dephasing_kraus;
	std::vector<Kraus> two_qubit_damping_kraus;
	std::vector<Kraus> two_qubit_dephasing_kraus;

	Circuit<qtraits> c;
	Circuit<qtraits> real_circuit;

	int one_qubit_gate_time;
	int two_qubit_gate_time;

	std::vector<int> T1;
	std::vector<int> T2;

	std::vector<fp_t> one_qubit_gate_error;
	std::vector<std::vector<fp_t>> two_qubit_gate_error;
	std::vector<std::vector<GateType>> clock_cycle;

	bool pre_gen_used = false;

	static struct Config {
		int one_qubit_gate_time = 1;
		int two_qubit_gate_time = 2;
		int global_T1 = 1000;
		int global_T2 = 1000;
		fp_t one_qubit_gate_error = 0.01;
		fp_t two_qubit_gate_error = 0.02;
	} default_config;

	RealCircuit() {}

	void assign_circuit(const Circuit<qtraits>& c_) {
		one_qubit_gate_error.assign(c_.max_qubit, default_config.one_qubit_gate_error);
		T1.assign(c_.max_qubit, default_config.global_T1);
		T2.assign(c_.max_qubit, default_config.global_T2);

		one_qubit_gate_time = default_config.one_qubit_gate_time;
		two_qubit_gate_time = default_config.two_qubit_gate_time;

		two_qubit_gate_error.resize(c_.max_qubit);
		for (auto& e : two_qubit_gate_error) {
			e.resize(c_.max_qubit, default_config.two_qubit_gate_error);
		}
		clock_cycle.resize(c_.max_qubit);
	}

	explicit RealCircuit(const Circuit<qtraits>& c_) :
		c(c_.valuecopy()) {

		assign_circuit(c_);
	}

	RealCircuit(const RealCircuit<qtraits>&) = delete;
	//RealCircuit(RealCircuit<qtraits>&& rc) {
	//	c = std::move(rc.c);
	//	real_circuit = std::move(rc.real_circuit);
	//	one_qubit_gate_time = rc.one_qubit_gate_time;
	//	two_qubit_gate_time = rc.two_qubit_gate_time;
	//	T1 = std::move(rc.T1);
	//	T2 = std::move(rc.T2);

	//	one_qubit_gate_error = std::move(rc.one_qubit_gate_error);
	//	two_qubit_gate_error = std::move(rc.two_qubit_gate_error);
	//	clock_cycle = std::move(rc.clock_cycle);
	//	pre_gen_used = rc.pre_gen_used;
	//}

	RealCircuit<qtraits>&& valuecopy() {
		RealCircuit<qtraits> rc;
		rc.c = c.valuecopy();
		rc.real_circuit = real_circuit.valuecopy();
		rc.one_qubit_gate_time = one_qubit_gate_time;
		rc.two_qubit_gate_time = two_qubit_gate_time;
		rc.T1 = T1;
		rc.T2 = T2;
		rc.one_qubit_gate_error = one_qubit_gate_error;
		rc.two_qubit_gate_error = two_qubit_gate_error;
		rc.clock_cycle = clock_cycle;
		rc.pre_gen_used = pre_gen_used;
		return rc;
	}

	void analyze_clock() {
		for (const auto &gate : c.gates) {
			assert(gate.controllers.size() == 0, "No controller is allowed.");
			assert(gate.dag == false, "Flatten first.");

			switch (gate.type) {
			case GateType::RX:
			case GateType::RY:
			case GateType::RZ:
			case GateType::H:
				add_single_qubit(gate, gate.qubits[0]);
				break;
			case GateType::CNOT:
				add_two_qubit(gate, gate.qubits[0], gate.qubits[1]);
				break;
			case GateType::I:
				add_I(gate, gate.qubits[0]);
				break;
			default:
				assert(false, "Not supported gate type.");
			}
		}
	}

	void generate_kraus_op() {

		qid qsize = c.max_qubit;
		/* one_qubit_damping_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t *kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_damping_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag,
				T1[i], one_qubit_gate_time);

			Kraus one_damping;
			one_damping.kraus0_real = kraus0_real;
			one_damping.kraus0_imag = kraus0_imag;
			one_damping.kraus1_real = kraus1_real;
			one_damping.kraus1_imag = kraus1_imag;

			one_qubit_damping_kraus.push_back(one_damping);
		}

		/* one_qubit_dephasing_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t *kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_dephasing_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag, T1[i],
				T2[i], one_qubit_gate_time);

			Kraus one_dephasing;
			one_dephasing.kraus0_real = kraus0_real;
			one_dephasing.kraus0_imag = kraus0_imag;
			one_dephasing.kraus1_real = kraus1_real;
			one_dephasing.kraus1_imag = kraus1_imag;

			one_qubit_dephasing_kraus.push_back(one_dephasing);
		}

		/* two_qubit_damping_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t *kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_damping_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag,
				T1[i], two_qubit_gate_time);

			Kraus two_damping;
			two_damping.kraus0_real = kraus0_real;
			two_damping.kraus0_imag = kraus0_imag;
			two_damping.kraus1_real = kraus1_real;
			two_damping.kraus1_imag = kraus1_imag;

			two_qubit_damping_kraus.push_back(two_damping);
		}

		/* two_qubit_dephasing_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t *kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t *kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_dephasing_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag, T1[i],
				T2[i], two_qubit_gate_time);

			Kraus two_dephasing;
			two_dephasing.kraus0_real = kraus0_real;
			two_dephasing.kraus0_imag = kraus0_imag;
			two_dephasing.kraus1_real = kraus1_real;
			two_dephasing.kraus1_imag = kraus1_imag;

			two_qubit_dephasing_kraus.push_back(two_dephasing);
		}
	}

	void add_single_qubit(Gate<qtraits> g, qid q) {
		clock_cycle[q].push_back(g.type);
		for (int i = 1; i < one_qubit_gate_time; ++i)
			clock_cycle[q].push_back(GateType(-1));
		real_circuit - g;
	}

	void add_two_qubit(Gate<qtraits> g, qid q1, qid q2) {
		size_t front1 = clock_cycle[q1].size();
		size_t front2 = clock_cycle[q2].size();

		if (front1 > front2) {
			clock_cycle[q2].push_back(GateType::I);
			for (int i = 1; i < front1 - front2; ++i)
				clock_cycle[q2].push_back(GateType(-1));
			real_circuit - I(q2, front1 - front2);
		}
		else if (front2 > front1) {
			clock_cycle[q1].push_back(GateType::I);
			for (int i = 1; i < front2 - front1; ++i)
				clock_cycle[q1].push_back(GateType(-1));
			real_circuit - I(q1, front2 - front1);
		}

		clock_cycle[q1].push_back(g.type);
		clock_cycle[q2].push_back(g.type);

		for (int i = 1; i < two_qubit_gate_time; ++i) {
			clock_cycle[q1].push_back(GateType(-1));
			clock_cycle[q2].push_back(GateType(-1));
		}
		real_circuit - g;
	}

	void add_I(Gate<qtraits> g, qid q) {
		clock_cycle[q].push_back(g.type);
		for (int i = 1; i < int(g.argument); ++i)
			clock_cycle[q].push_back(GateType(-1));
		real_circuit - g;
	}

	void ready() {
		analyze_clock();
		if (use_pre_gen && !pre_gen_used) {
			pre_gen_used = true;
			one_qubit_damping_kraus.assign(pre_gen_1qdam.begin(), pre_gen_1qdam.end());
			one_qubit_dephasing_kraus.assign(pre_gen_1qdep.begin(), pre_gen_1qdep.end());
			two_qubit_damping_kraus.assign(pre_gen_2qdam.begin(), pre_gen_2qdam.end());
			two_qubit_dephasing_kraus.assign(pre_gen_2qdep.begin(), pre_gen_2qdep.end());
		}
		else {
			generate_kraus_op();
		}
	}

	std::string to_string() const {
		std::stringstream ss;
		for (auto cc : clock_cycle) {
			for (auto p : cc) {
				if ((int)p == -1)
					ss << std::setw(3) << " ";
				else
					ss << std::setw(3) << (int)p;
			}
			ss << std::endl;
		}
		return ss.str();
	}

	~RealCircuit() {
		if (pre_gen_used)
			return;
		for (auto p : one_qubit_damping_kraus) {
			free(p.kraus0_real);
			free(p.kraus0_imag);
			free(p.kraus1_real);
			free(p.kraus1_imag);
		}
		for (auto p : one_qubit_dephasing_kraus) {
			free(p.kraus0_real);
			free(p.kraus0_imag);
			free(p.kraus1_real);
			free(p.kraus1_imag);
		}
		for (auto p : two_qubit_damping_kraus) {
			free(p.kraus0_real);
			free(p.kraus0_imag);
			free(p.kraus1_real);
			free(p.kraus1_imag);
		}
		for (auto p : two_qubit_dephasing_kraus) {
			free(p.kraus0_real);
			free(p.kraus0_imag);
			free(p.kraus1_real);
			free(p.kraus1_imag);
		}
	}

	static std::vector<Kraus> pre_gen_1qdam;
	static std::vector<Kraus> pre_gen_1qdep;
	static std::vector<Kraus> pre_gen_2qdam;
	static std::vector<Kraus> pre_gen_2qdep;

	static void pre_gen_kraus(size_t qsize) {
		use_pre_gen = true;
		auto& T1 = decltype(*this)::default_config.global_T1;
		auto& T2 = decltype(*this)::default_config.global_T2;
		auto& one_qubit_gate_time = decltype(*this)::default_config.one_qubit_gate_time;
		auto& two_qubit_gate_time = decltype(*this)::default_config.two_qubit_gate_time;

		/* one_qubit_damping_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t* kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_damping_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag,
				T1, one_qubit_gate_time);

			Kraus one_damping;
			one_damping.kraus0_real = kraus0_real;
			one_damping.kraus0_imag = kraus0_imag;
			one_damping.kraus1_real = kraus1_real;
			one_damping.kraus1_imag = kraus1_imag;

			pre_gen_1qdam.push_back(one_damping);
		}

		/* one_qubit_dephasing_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t* kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_dephasing_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag, T1,
				T2, one_qubit_gate_time);

			Kraus one_dephasing;
			one_dephasing.kraus0_real = kraus0_real;
			one_dephasing.kraus0_imag = kraus0_imag;
			one_dephasing.kraus1_real = kraus1_real;
			one_dephasing.kraus1_imag = kraus1_imag;

			pre_gen_1qdep.push_back(one_dephasing);
		}

		/* two_qubit_damping_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t* kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_damping_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag,
				T1, two_qubit_gate_time);

			Kraus two_damping;
			two_damping.kraus0_real = kraus0_real;
			two_damping.kraus0_imag = kraus0_imag;
			two_damping.kraus1_real = kraus1_real;
			two_damping.kraus1_imag = kraus1_imag;

			pre_gen_2qdam.push_back(two_damping);
		}

		/* two_qubit_dephasing_kraus */
		for (qid i = 0; i < qsize; ++i) {
			fp_t* kraus0_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_real = (fp_t*)malloc(4 * sizeof(fp_t));
			fp_t* kraus1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

			get_dephasing_kraus_op(kraus0_real, kraus0_imag,
				kraus1_real, kraus1_imag, T1,
				T2, two_qubit_gate_time);

			Kraus two_dephasing;
			two_dephasing.kraus0_real = kraus0_real;
			two_dephasing.kraus0_imag = kraus0_imag;
			two_dephasing.kraus1_real = kraus1_real;
			two_dephasing.kraus1_imag = kraus1_imag;

			pre_gen_2qdep.push_back(two_dephasing);
		}
	}
	static void set_use_pre_gen(bool use) {
		use_pre_gen = use;
	}
	static bool use_pre_gen;

};

template<typename qtraits>
bool RealCircuit<qtraits>::use_pre_gen = false;

template<typename qtraits>
typename RealCircuit<qtraits>::Config 
RealCircuit<qtraits>::default_config;

template<typename qtraits>
std::vector<typename RealCircuit<qtraits>::Kraus> 
RealCircuit<qtraits>::pre_gen_1qdam;

template<typename qtraits>
std::vector<typename RealCircuit<qtraits>::Kraus>
RealCircuit<qtraits>::pre_gen_1qdep;

template<typename qtraits>
std::vector<typename RealCircuit<qtraits>::Kraus>
RealCircuit<qtraits>::pre_gen_2qdam;

template<typename qtraits>
std::vector<typename RealCircuit<qtraits>::Kraus>
RealCircuit<qtraits>::pre_gen_2qdep;

template<typename qtraits_t = default_qtraits>
struct result_analyzer {
	using qtraits = qtraits_t;
	using qid = typename qtraits::qidx_t;
	using fp_t = typename qtraits::value_t;
	using uint_t = typename qtraits::idx_t;

	static size_t count_zero(uint_t* result, size_t shots) {
		size_t zeros = 0;
		for (size_t i = 0; i < shots; ++i) {
			if (result[i] == 0)
				zeros++;
		}
		return zeros;
	}

	static size_t count_n(uint_t* result, size_t shots, uint_t n) {
		size_t count = 0;
		for (size_t i = 0; i < shots; ++i) {
			if (result[i] == n)
				count++;
		}
		return count;
	}

	static std::vector<qid> generate_qubit_idxes(qid max_qubit) {
		std::vector<qid> q(max_qubit, 0);
		for (qid i = 0; i < max_qubit; ++i) {
			q[i] = i;
		}
		return q;
	}

	static std::pair<uint_t, size_t> find_max(uint_t* result, size_t shots) {
		std::map<uint_t, size_t> max_counter;
		std::pair<uint_t, size_t> max = { 0,0 };
		for (size_t i = 0; i < shots; ++i) {
			if (max_counter.find(result[i]) == max_counter.end()) {
				max_counter[result[i]] = 0;
			}
			max_counter[result[i]]++;
			size_t &p = max_counter[result[i]];
			if (p > max.second) {
				max = { result[i], p };
			}
		}
		return max;
	}

	static fp_t get_expectation_from_amplitude(fp_t *real, fp_t *imag, size_t size, fp_t *diag_h) {
		fp_t expectation = 0;
		for (size_t i = 0; i < size; ++i) {
			fp_t prob = real[i] * real[i] + imag[i] * imag[i];
			expectation += prob * diag_h[i];
		}
		return expectation;
	}

	static double get_expectation(const uint_t *res, const uint_t size, const std::vector<double> &diag) {
		double exp = 0;
		for (uint_t i = 0; i < size; ++i) {
			exp += diag[res[i]];
		}
		return exp / size;
	}

	static std::vector<fp_t> get_meas_probs(const uint_t *res, const uint_t shots, const uint_t size) {
		std::vector<fp_t> probs(size, 0);
		for (uint_t i = 0; i < shots; ++i) {
			probs[res[i]] += (1.0 / shots);
		}
		return probs;
	}

	static std::vector<uint_t> get_meas_count(const uint_t *res, const uint_t shots, const uint_t size) {
		std::vector<uint_t> count(size, 0);
		for (uint_t i = 0; i < shots; ++i) {
			count[res[i]]++;
		}
		return count;
	}

	static bool check_error(fp_t *prob1, fp_t *prob2, uint_t size, fp_t error_bound) {
		for (uint_t i = 0; i < size; ++i) {
			if (abs(prob1[i] - prob2[i]) > error_bound) {
				return false;
			}
		}
		return true;
	}

	static fp_t get_norm_2(const fp_t *prob1, const fp_t *prob2, const uint_t size) {
		fp_t norm2 = 0;
		for (uint_t i = 0; i < size; ++i) {
			norm2 += (prob1[i] - prob2[i])*(prob1[i] - prob2[i]);
		}
		norm2 = sqrt(norm2);
		return norm2;
	}

	static fp_t get_norm_inf(const fp_t *prob1, const fp_t *prob2, const uint_t size) {
		fp_t norminf = 0;
		for (uint_t i = 0; i < size; ++i) {
			fp_t ninf = abs(prob1[i] - prob2[i]);
			if (ninf > norminf) norminf = ninf;
		}
		return norminf;
	}
};

template<typename qtraits = default_qtraits>
struct simulator_v1 {
	using fp_t = typename qtraits::value_t;
	using uint_t = typename qtraits::idx_t;
	using qid = typename qtraits::qidx_t;

	static void run_circuit(const Circuit<qtraits> &c, fp_t *real, fp_t *imag, uint_t size) {
		const qid &qn = c.max_qubit;
		// auto &groupdefs = c.groupdefs;
		size_t groupdef_iter = 0;

		for (const Gate<qtraits> &gate : c.gates) {
			switch (gate.type) {
			case GateType::RX:
				state_manipulator<qtraits>::rx(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::RY:
				state_manipulator<qtraits>::ry(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::RZ:
				state_manipulator<qtraits>::rz(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::H:
				state_manipulator<qtraits>::h(real, imag, size, gate.qubits[0]);
				break;
			case GateType::CNOT:
				state_manipulator<qtraits>::cnot(real, imag, size, gate.qubits[0], gate.qubits[1]);
				break;
			case GateType::I:
				break;
			default:
				assert(false, "Bad Type.");
			}
		}
	}

	static uint_t *simulate_N(size_t n, const Circuit<qtraits> &c, RandomEngine *rng) {
		// first allocate n times memory
		uint_t size = pow2(c.max_qubit);
		fp_t *real_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		fp_t *imag_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		uint_t *result = (uint_t*)malloc(n * sizeof(uint_t));
		memset(result, 0, sizeof(uint_t)*n);

		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;

			run_circuit(c, real, imag, size);

			fp_t randnum = (fp_t)(*rng)();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
		free(real_n);
		free(imag_n);
		return result;
	}

	/* 
		This version will assume the real, imag is refreshed to zero state. 
		and result has enough space. 
	*/
	static void simulate_N(size_t n, const Circuit<qtraits> &c, RandomEngine** rng, 
		fp_t *real_n, fp_t *imag_n, uint_t *result) {

		uint_t size = pow2(c.max_qubit);
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;

			run_circuit(c, real, imag, size);

			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
	}

	static uint_t* simulate_N_threads(size_t n, const Circuit<qtraits>& c, RandomEngine** rng) {
		// first allocate n times memory
		uint_t size = pow2(c.max_qubit);
		fp_t* real_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		fp_t* imag_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		uint_t* result = (uint_t*)malloc(n * sizeof(uint_t));
		memset(result, 0, sizeof(uint_t) * n);

#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];

			real[0] = 1;

			run_circuit(c, real, imag, size);

			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
		free(real_n);
		free(imag_n);
		return result;
	}

	static void simulate_N_threads(size_t n, const Circuit<qtraits> &c, RandomEngine** rng,
		fp_t* real_n, fp_t* imag_n, uint_t* result) {

		uint_t size = pow2(c.max_qubit);
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];

			real[0] = 1;	
			run_circuit(c, real, imag, size);

			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
	}

	static void run_real_circuit(const RealCircuit<qtraits> &rc,
		fp_t *real, fp_t *imag, uint_t size, RandomEngine* rng) {

		auto &p1 = rc.one_qubit_gate_error;
		auto &p2 = rc.two_qubit_gate_error;

		auto &T1 = rc.T1;
		auto &T2 = rc.T2;
		auto &one_qubit = rc.one_qubit_gate_time;
		auto &two_qubit = rc.two_qubit_gate_time;

		const Circuit<qtraits>& c = rc.real_circuit;
		const qid &qn = c.max_qubit;

		/*debug_info_s(Circuit:);
		debug_print_circuit;*/

		for (auto &gate : c.gates) {
			// incoherent error
			// for first qubit
			fp_t *damp0_real = nullptr, *damp0_imag = nullptr,
				*damp1_real = nullptr, *damp1_imag = nullptr;
			fp_t *dephase0_real = nullptr, *dephase0_imag = nullptr,
				*dephase1_real = nullptr, *dephase1_imag = nullptr;
			// for second qubit
			fp_t *damp0_real2 = nullptr, *damp0_imag2 = nullptr,
				*damp1_real2 = nullptr, *damp1_imag2 = nullptr;
			fp_t *dephase0_real2 = nullptr, *dephase0_imag2 = nullptr,
				*dephase1_real2 = nullptr, *dephase1_imag2 = nullptr;

			if (gate.type == GateType::I) {
				damp0_real = (fp_t*)malloc(4 * sizeof(fp_t));
				damp0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
				damp1_real = (fp_t*)malloc(4 * sizeof(fp_t));
				damp1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

				dephase0_real = (fp_t*)malloc(4 * sizeof(fp_t));
				dephase0_imag = (fp_t*)malloc(4 * sizeof(fp_t));
				dephase1_real = (fp_t*)malloc(4 * sizeof(fp_t));
				dephase1_imag = (fp_t*)malloc(4 * sizeof(fp_t));

				get_damping_kraus_op(
					damp0_real, damp0_imag, damp1_real, damp1_imag,
					T1[gate.qubits[0]], (int)gate.argument);

				get_dephasing_kraus_op(
					dephase0_real, dephase0_imag, dephase1_real, dephase1_imag,
					T1[gate.qubits[0]], T2[gate.qubits[0]], (int)gate.argument);
			}
			else if (gate.qubits.size() == 1) {
				auto &damp = rc.one_qubit_damping_kraus[gate.qubits[0]];
				damp0_real = damp.kraus0_real;
				damp0_imag = damp.kraus0_imag;
				damp1_real = damp.kraus1_real;
				damp1_imag = damp.kraus1_imag;

				auto &dephase = rc.one_qubit_dephasing_kraus[gate.qubits[0]];
				dephase0_real = dephase.kraus0_real;
				dephase0_imag = dephase.kraus0_imag;
				dephase1_real = dephase.kraus1_real;
				dephase1_imag = dephase.kraus1_imag;
			}
			else if (gate.qubits.size() == 2) {
				auto &damp = rc.two_qubit_damping_kraus[gate.qubits[0]];
				damp0_real = damp.kraus0_real;
				damp0_imag = damp.kraus0_imag;
				damp1_real = damp.kraus1_real;
				damp1_imag = damp.kraus1_imag;

				auto &dephase = rc.two_qubit_dephasing_kraus[gate.qubits[0]];
				dephase0_real = dephase.kraus0_real;
				dephase0_imag = dephase.kraus0_imag;
				dephase1_real = dephase.kraus1_real;
				dephase1_imag = dephase.kraus1_imag;

				auto &damp2 = rc.two_qubit_damping_kraus[gate.qubits[1]];
				damp0_real2 = damp.kraus0_real;
				damp0_imag2 = damp.kraus0_imag;
				damp1_real2 = damp.kraus1_real;
				damp1_imag2 = damp.kraus1_imag;

				auto &dephase2 = rc.two_qubit_dephasing_kraus[gate.qubits[1]];
				dephase0_real2 = dephase.kraus0_real;
				dephase0_imag2 = dephase.kraus0_imag;
				dephase1_real2 = dephase.kraus1_real;
				dephase1_imag2 = dephase.kraus1_imag;
			}
			else assert(false, "Bad Gate Number.");
			double r0 = (*rng)();
			double r1 = (*rng)();

		/*	debug_info_s(Before Kraus);
			debug_output_state;
			debug_display(r0);
			debug_display(r1);*/

			state_manipulator<qtraits>::perform_kraus(real, imag, size, damp0_real, damp0_imag,
				damp1_real, damp1_imag, gate.qubits[0], r0);

			state_manipulator<qtraits>::perform_kraus(real, imag, size, dephase0_real, dephase0_imag,
				dephase1_real, dephase1_imag, gate.qubits[0], r1);

			if (gate.type == GateType::I) {
				free(damp0_real);
				free(damp0_imag);
				free(damp1_real);
				free(damp1_imag);
				free(dephase0_real);
				free(dephase0_imag);
				free(dephase1_real);
				free(dephase1_imag);
			}

			if (gate.qubits.size() == 2) {
				double p0 = (*rng)();
				double p1 = (*rng)();
				state_manipulator<qtraits>::perform_kraus(real, imag, size, damp0_real2, damp0_imag2,
					damp1_real2, damp1_imag2, gate.qubits[1], p0);
				state_manipulator<qtraits>::perform_kraus(real, imag, size, dephase0_real2, dephase0_imag2,
					dephase1_real2, dephase1_imag2, gate.qubits[1], p1);
			}

			// coherent error
			double p_bad_gate = (*rng)();
			if (gate.type != GateType::I) {
				if (gate.qubits.size() == 1) {
					if (p_bad_gate < p1[gate.qubits[0]]) {
						goto SkipGate;
					}
				}
				else if (gate.qubits.size() == 2) {
					if (p_bad_gate < p2[gate.qubits[0]][gate.qubits[1]]) {
						goto SkipGate;
					}
				}
				else assert(false, "Bad Gate Number.");
			}

			/*debug_info_s(After Kraus);
			debug_output_state;*/

			switch (gate.type) {
			case GateType::RX:
				state_manipulator<qtraits>::rx(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::RY:
				state_manipulator<qtraits>::ry(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::RZ:
				state_manipulator<qtraits>::rz(real, imag, size, gate.argument, gate.qubits[0]);
				break;
			case GateType::H:
				state_manipulator<qtraits>::h(real, imag, size, gate.qubits[0]);
				break;
			case GateType::CNOT:
				state_manipulator<qtraits>::cnot(real, imag, size, gate.qubits[0], gate.qubits[1]);
				break;
			case GateType::I:
				break;
			default:
				assert(false, "Bad Type.");
			}
		SkipGate:;
			/*debug_info_s(Last);
			debug_output_state;
			debug_pause;*/
		}
	}

	static uint_t *simulate_N_noisy(size_t n, const RealCircuit<qtraits> &c, RandomEngine **rng) {
		// first allocate n times memory
		uint_t size = pow2(c.real_circuit.max_qubit);

		fp_t *real_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		fp_t *imag_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		uint_t *result = (uint_t*)malloc(n * sizeof(uint_t));
		memset(result, 0, sizeof(uint_t)*n);

		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;
			run_real_circuit(c, real, imag, size, rng[i]);
			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}

		free(real_n);
		free(imag_n);

		return result;
	}

	static void simulate_N_noisy(size_t n, const RealCircuit<qtraits>& c, RandomEngine** rng,
		fp_t *real_n, fp_t *imag_n, uint_t* result) {
		// first allocate n times memory
		uint_t size = pow2(c.real_circuit.max_qubit);

		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;
			run_real_circuit(c, real, imag, size, rng[i]);
			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
	}

	static uint_t *simulate_N_threads_noisy(size_t n, const RealCircuit<qtraits> &c, RandomEngine **rng) {
		// first allocate n times memory
		uint_t size = pow2(c.real_circuit.max_qubit);
		fp_t *real_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		fp_t *imag_n = (fp_t*)malloc(n * size * sizeof(fp_t));
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

		uint_t *result = (uint_t*)malloc(n * sizeof(uint_t));
		// memset(result, 0, sizeof(uint_t)*n);

#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;

			run_real_circuit(c, real, imag, size, rng[i]);
			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
		free(real_n);
		free(imag_n);

		return result;
	}

	static void simulate_N_threads_noisy(size_t n, const RealCircuit<qtraits>& c, RandomEngine** rng,
		fp_t* real_n, fp_t* imag_n, uint_t *result) {
		// first allocate n times memory
		uint_t size = pow2(c.real_circuit.max_qubit);
		memset(real_n, 0, n * size * sizeof(fp_t));
		memset(imag_n, 0, n * size * sizeof(fp_t));

#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			fp_t* real = &real_n[i * size];
			fp_t* imag = &imag_n[i * size];
			real[0] = 1;

			run_real_circuit(c, real, imag, size, rng[i]);
			fp_t randnum = (fp_t)(*rng[i])();
			fp_t total_prob = 0;
			for (uint_t j = 0; j < size; ++j) {
				fp_t pi = real[j] * real[j] + imag[j] * imag[j];
				total_prob += pi;
				if (randnum < total_prob) {
					result[i] = j; break;
				}
			}
		}
	}
};

/* simulator modes */
constexpr bool noisy = true;
constexpr bool noise_free = false;
constexpr bool single_thread = false;
constexpr bool multi_threads = true;

ns_end