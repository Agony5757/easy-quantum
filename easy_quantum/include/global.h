#pragma once

#include <vector>
#include <algorithm>
#include <type_traits>
#include <string>
#include <sstream>
#include <array>
#include <functional>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include <map>
#include <unordered_map>
#include <any>
#include <optional>
#include <mutex>
#include <complex>

#define ns_easyquantum namespace easyquantum {
#define ns_end	 }
#define using_ns_easyquantum using namespace easyquantum;

template<typename value_type, typename idx_type, typename qidx_type>
struct quantum_traits {
	using value_t = value_type;
	using idx_t = idx_type;
	using qidx_t = qidx_type;
};

using default_qtraits = quantum_traits<double, uint64_t, int>;

/* For loops */
#define For_Size(container) container.size()
#define For(container, elem, start, stop) for(size_t elem = start; elem < stop; ++ elem)
#define ForAll(container, elem) For(container, elem, 0, For_Size(container))
#define ForAlli(container) ForAll(container,i)
#define Loop(N, i) for (size_t i=0;i<N;++i)
#define LoopA(N) Loop(N, ____slslslsls____)
#define LoopI(N) Loop(N, i)