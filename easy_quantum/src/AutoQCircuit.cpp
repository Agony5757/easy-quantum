#include <AutoQCircuit.h>

ns_easyquantum

std::map<size_t, qubitpool> qubitpool::_qubitpools;
size_t qubitpool::context_id = 0;

ns_end