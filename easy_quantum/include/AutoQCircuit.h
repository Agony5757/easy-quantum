#pragma once
#include <global.h>
#include <vector>
#include <exception>

ns_easyquantum

struct QubitStatus {
	int id;
	bool allocated;
};

struct CMemStatus {
	int id;
	bool allocated;
	int value;
};

class qsystem {
private:
	std::vector<QubitStatus> qubits;
	std::vector<CMemStatus> CMemStatus;

	static qsystem system;
public:
	static qsystem& instance() { return system; }
	int init_system(int nqubits);
	int first_free_qubit();
	int first_free_cmem();
	bool allocate_qubit(int id);
	bool allocate_cbit(int id);
};

class qubit {
	// construct to allocate
	// destruct to free
	int id;
	qubit() {
		
		qsystem& sys = qsystem::instance();
		int id = sys.first_free_qubit();
		if (id >= 0) {
			sys.allocate_qubit(id);
		}
		else {
			throw std::runtime_error("No more free qubits");
		}

	}

	qubit(const qubit&) = delete;
};

class qubits {

};

class cbit {
	cbit() {

	}
};

class qmodule {
	std::vector<std::vector<int>> n_wires;
};

ns_end // ns_easyquantum
