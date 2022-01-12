#include <global.h>
#include <mathutil.h>
#include <Windows.h>

ns_easyquantum

template<typename Ty>
struct safe_allocator {
	static_assert(std::is_floating_point_v<Ty> && !std::is_const_v<Ty> && !std::is_volatile_v<Ty>,
		"Class safe_allocator only supports floating points value types.");

	static bool has_checked_memory;
	static size_t memory_length;
	static size_t max_qubit_number;

	static void memory_check() {
		MEMORYSTATUS ms;
		::GlobalMemoryStatus(&ms);
		
		// get available memory count, in bytes
		SIZE_T phymemory = ms.dwAvailPhys;

		// Ty = double 
		// 10 qubit 16KB
		// 20 qubit 16MB
		// 30 qubit 16GB
		// 40 qubit 16TB
		memory_length = (size_t)phymemory;
		has_checked_memory = true;

		max_qubit_number = log2(phymemory / 10);		
	}

	static size_t get_max_qubit_number() {
		if (!has_checked_memory)
			memory_check();

		return max_qubit_number;
	}

	Ty* allocate(const size_t _Count) {
		if (!has_checked_memory)
			memory_check();

		assert(_Count <= pow2(max_qubit_number) * 10,
			"Can only allocate 10*2^max_qubit numbers at once.");

		return new Ty[_Count];
	}

	void deallocate(Ty* const _Ptr, const size_t _Count) {
		delete[] _Ptr;
	}
};

template<typename Ty>
bool safe_allocator<Ty>::has_checked_memory = false;
template<typename Ty>
size_t safe_allocator<Ty>::memory_length = 1;
template<typename Ty>
size_t safe_allocator<Ty>::max_qubit_number = 1;

ns_end