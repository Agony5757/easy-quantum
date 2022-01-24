#include <iostream>
#include <easy_quantum>

using_ns_easyquantum

using namespace std;

void test1() {
	using qtraits = default_qtraits;
		
	Circuit c;
	c - H(1);
	c - H(2);
	c - H(3);
	c - H(4);
	c - H(0);
	RealCircuit rc(c);
	rc.ready();

	constexpr int shots = 100;

	DefaultRandomEngine rng[1];

	auto result = simulator_v1<qtraits>::simulate_N(shots, c, rng);

	RandomEngine** rngs = make_engines(shots);
	auto result2 = simulator_v1<qtraits>::simulate_N_noisy(shots, rc, rngs);
	free_engines(rngs, shots);
	
	LoopI(shots) {
		cout << i << " : "<<result[i] << endl;
	}
}


void test4() {

	using qid = typename default_qtraits::qidx_t;
	Hamiltonian<default_qtraits> h1(3.5, { { 'x', 1 } });
	Hamiltonian<default_qtraits> h2({ { 'x',1 }, {'y',2} });
	Hamiltonian<default_qtraits> h3({ 3.5, { { 'x',1 }, {'y',2} } });
	Hamiltonian<default_qtraits> h4(
	{
		{ { { 'x',1 }, {'y',2} }, 3.5},
		{ { { 'x',1 }, {'y',3} }, 3.5} 
	});

	cout << h1.to_string() << endl;
	cout << h2.to_string() << endl;
	cout << h3.to_string() << endl;
	cout << h4.to_string() << endl;

}

void test5() {
	using qid = typename default_qtraits::qidx_t;
	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',1 }, {'z',2} }, 3.5},
			{ { { 'z',1 }, {'z',3} }, 3.5}
		});

	bool suc = h.get_full_diagonal({ 1,2,3 });
	if (suc) {
		cout << h.full_diagonal.to_string();
	}
	else {
		cout << "Fail";
	}
}

int main() {
	TestSet::getInstance().run_all_tests();
	getchar();
	return 0;
}