#include <easy_quantum>

using_ns_easyquantum

Test(HamiltonianTest1) {
	TestResult res;
	using qid = typename default_qtraits::qidx_t;
	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',1 }, {'z',2} }, 3.5},
			{ { { 'z',1 }, {'z',3} }, 3.5}
		});

	using vector_t = decltype(h)::vector_t;

	bool suc = h.get_full_diagonal({ 1,2,3 });
	res &= testcheck(suc, "Fail getting full diagonal.");
	vector_t& diag = h.full_diagonal;
	vector_t checkdiag = { 7,-7,0,0,0,0,-7,7 };
	res &= testcheck(diag == checkdiag, "Bad value");

	return res;
}

