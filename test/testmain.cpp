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

template<typename qtraits>
Circuit<qtraits> getcircuit(Variable<qtraits> &a) {
	Circuit<qtraits> cir;
	cir - RX(0, a)
		- RX(1, a)
		- RX(2, a)
		- RX(3, 0)
		- RX(4, a);

	return cir;
}

void test3() {
	Variable<default_qtraits> a(20);	
	Circuit cir = getcircuit<default_qtraits>(a);
	Circuit cir2;

	cir2 - cir.valuecopy()
		- RX(1, a);
	cout << cir.to_qasm();
	cout << cir2.to_qasm();
	a.set_value(50);
	cout << cir.to_qasm();
	cout << cir2.to_qasm();
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

void test6() {
	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;
	Hamiltonian<default_qtraits> h(
		{
			{ { { 'x',1 }, {'z',2} }, 3.5},
			{ { { 'z',1 }, {'z',3} }, 3.5}
		});

	//c.max_qubit = 4;
	Circuit<default_qtraits> c;
	c - RX(0, 2.3);
	c - RX(1, 4.0);
	c - RX(2, 1.7);
	c - RX(3, 0.8);
	ObjectiveFunction obj({}, c, h);

	fp_t exp1 = obj.get_expectation_montecarlo(100000, single_thread, noise_free);
	cout << "expectation (monte) (st) = " << exp1 << endl;

	fp_t exp2 = obj.get_expectation_montecarlo(100000, multi_threads, noise_free);
	cout << "expectation (monte) (mt) = " << exp2 << endl;

	fp_t exp3 = obj.get_expectation_prob_noisefree();
	cout << "expectation (probs) (st) = " << exp3 << endl;
}

void test7() {
	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;
	Hamiltonian<default_qtraits> h(
		{
			{ { { 'x',1 }, {'z',2} }, 3.5},
			{ { { 'z',1 }, {'z',3} }, 3.5}
		});
	Variable<default_qtraits> v;
	Circuit<default_qtraits> c;
	c - RX(0, v);
	c - RX(1, 4.0);
	c - RX(2, 1.7);
	c - RX(3, 0.8);
	//c.max_qubit = 4;

	ObjectiveFunction obj({v}, c, h);
	v.set_value(2.3);

	fp_t exp = obj.get_expectation_prob_noisefree();
	cout << "expectation (probs) = " << exp << endl;

	vector<fp_t> der = obj.get_derivatives_prob_noisefree();
	cout << "derivative (probs) = " << der[0] << endl;
}

void test7_1() {
	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;
	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',0 }, }, 15},
		});
	Variable<default_qtraits> v;
	Circuit<default_qtraits> c;
	c - RX(0, v);
	//c.max_qubit = 4;

	ObjectiveFunction obj({ v }, c, h);
	v.set_value(pi);

	fp_t exp = obj.get_expectation_prob_noisefree();
	cout << "expectation (probs) = " << exp << endl;

	vector<fp_t> der = obj.get_derivatives_prob_noisefree();
	cout << "derivative (probs) = " << der[0] << endl;
}

int main() {
	TestSet::getInstance().run_all_tests();
	getchar();
	return 0;
}