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
};

template<size_t pool_id>
struct pool {
	static std::map<size_t, pool<pool_id>> _pools;
	static size_t _context_id;		
	size_t pool_id;
	size_t max_num;
	std::vector<size_t> allocated;

	// can only called by the system::init()
	static void init_pool_system() {
		switch_context(0);
	}

	static void switch_context(size_t ctx, size_t max_num = 30) 
	{
		auto iter = _pools.find(ctx);
		if (iter == _pools.end()) {
			// create new pool
			pool& newpool = _pools[ctx] = pool();
			newpool.pool_id = ctx;
			newpool.max_num = max_num;
			newpool.allocated.resize(0, max_num);
			_context_id = ctx;
		}
		else {
			_context_id = ctx;
		}
	}

	size_t allocate() 
	{
		for (size_t i = 0; i < allocated.size(); ++i) {
			if (allocated[i] == 0) { return i; }
		}
		throw std::runtime_error("No free space.");
	}

	void copy_construct(size_t id)
	{
		allocated[id] += 1;
	}

	void destruct(size_t id)
	{
		allocated[id] += 1;
	}

	static pool& instance() {
		return _pools[_context_id];
	}

	static pool& instance(size_t pid) {
		return _pools[pid];
	}
};

template<size_t i>
static std::map<size_t, pool<i>> pool<i>::_pools;

template<size_t i>
static size_t pool<i>::_context_id = 0;

using qubitpool = pool<1>;
using cbitpool = pool<2>;

class qubit {
	// construct to allocate
	// destruct to free
	size_t id;
	size_t pool_id;
	qubit() {		
		qubitpool& pool = qubitpool::instance();
		id = pool.allocate();
		pool_id = pool.pool_id;
	}

	qubit(const qubit& old) {
		id = old.id;
		pool_id = old.id;
		qubitpool::instance(pool_id).copy_construct(id);
	}

	~qubit() {
		qubitpool::instance(pool_id).destruct(id);
	}

	qubit(const qubit&) = delete;
};

class cbit {
	cbit() {

	}
};

class qmodule {
	size_t pool_ctx_id;
	
};

ns_end // ns_easyquantum
