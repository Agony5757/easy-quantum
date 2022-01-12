#include <easy_quantum>

using_ns_easyquantum

Test(HamiltonianAndOperatorTest1) {
	TestResult res;
	
	Operator<default_qtraits> op1({ { 1,'z' } });
	Operator<default_qtraits> op2({ { { 1,'z' } }, 44 });
	Operator<default_qtraits> op3 = op1 + op2;
	res &= testcheck(op3.terms[{ { 1, 'z' } }] == 45., "Step 1");
	res &= testcheck(op3.terms[{ { 2, 'z' } }] == 0., "Step 2");

	return res;
}

Test(DifferentialTest1) {
	TestResult res;

	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;

	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',0 }, }, 15},
		});

	Variable<default_qtraits> v;
	Circuit<default_qtraits> c;
	c - RX(0, v);
	
	ObjectiveFunction obj({ v }, c, h);
	fp_t exp = -999, dif = -999;
	/*  H = 15 Z0
		|s(t)> = Rx(t)|0>

		E(pi/2) = 0;	E'(pi/2) = -15
		E(-pi/2) = 0;	E'(-pi/2) = 15
		E(pi) = -15;	E'(pi) = 0
		E(0) = 15;		E'(0) = 0;
	
	*/
	v.set_value(0);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 15., quantum_default_threshold), "Step 1/4");
	res &= testcheck(roundn(dif, 0., quantum_default_threshold), "Step 1/4");

	v.set_value(pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 0., quantum_default_threshold), "Step 2/4");
	res &= testcheck(roundn(dif, -15., quantum_default_threshold), "Step 2/4");

	v.set_value(-pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 0., quantum_default_threshold), "Step 3/4");
	res &= testcheck(roundn(dif, 15., quantum_default_threshold), "Step 3/4");

	v.set_value(pi);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, -15., quantum_default_threshold), "Step 4/4");
	res &= testcheck(roundn(dif, 0., quantum_default_threshold), "Step 4/4");

	return res;
}

Test(DifferentialTest2) {
	TestResult res;

	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;

	Hamiltonian<default_qtraits> h(
		{
			{ { { 'y',0 }, }, 15},
		});

	Variable<default_qtraits> v;
	Circuit<default_qtraits> c;
	c - RX(0, v);

	ObjectiveFunction obj({ v }, c, h);
	fp_t exp = -999, dif = -999;
	/*	H = 15 Y0
		|s(t)> = Rx(t)|0>

		E(0) = 0;		E'(0) = -15;
		E(pi/2) = -15;	E'(pi/2) = 0;
		E(-pi/2) = 15;	E'(-pi/2) = 0;
		E(pi) = 0;		E'(pi) = 15;		
	*/
	v.set_value(0);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 0., quantum_default_threshold), "Step 1/4");
	res &= testcheck(roundn(dif, -15., quantum_default_threshold), "Step 1/4");

	v.set_value(pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, -15., quantum_default_threshold), "Step 2/4");
	res &= testcheck(roundn(dif, 0., quantum_default_threshold), "Step 2/4");

	v.set_value(-pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 15., quantum_default_threshold), "Step 3/4");
	res &= testcheck(roundn(dif, 0., quantum_default_threshold), "Step 3/4");

	v.set_value(pi);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 0., quantum_default_threshold), "Step 4/4");
	res &= testcheck(roundn(dif, 15., quantum_default_threshold), "Step 4/4");

	return res;
}

Test(DifferentialTest3) {
	TestResult res;

	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;

	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',0 }, }, 15},
			{ { { 'y',0 }, }, 15},
		});

	Variable<default_qtraits> v;
	Circuit<default_qtraits> c;
	c - RX(0, v);
	ObjectiveFunction obj({ v }, c, h);

	fp_t exp = -999, dif = -999;
	/*
		H = 15 (Y0 + Z0)
		|s(t)> = Rx(t)|0>

		E(0) = 15;		E'(0) = -15;
		E(pi/2) = -15;	E'(pi/2) = -15
		E(-pi/2) = 15;	E'(-pi/2) = 15
		E(pi) = -15;	E'(pi) = 15

	*/
	v.set_value(0);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 15., quantum_default_threshold), "Step 1/4");
	res &= testcheck(roundn(dif, -15., quantum_default_threshold), "Step 1/4");

	v.set_value(pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, -15., quantum_default_threshold), "Step 2/4");
	res &= testcheck(roundn(dif, -15., quantum_default_threshold), "Step 2/4");

	v.set_value(-pi / 2);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, 15., quantum_default_threshold), "Step 3/4");
	res &= testcheck(roundn(dif, 15., quantum_default_threshold), "Step 3/4");

	v.set_value(pi);
	exp = obj.get_expectation_prob_noisefree();
	dif = obj.get_derivatives_prob_noisefree()[0];

	res &= testcheck(roundn(exp, -15., quantum_default_threshold), "Step 4/4");
	res &= testcheck(roundn(dif, 15., quantum_default_threshold), "Step 4/4");

	return res;
}

Test(DifferentialTest4) {
	TestResult res;

	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;
	
	/* H = 15 z0z1z3 + 78 x1y2z3 - 95 z1x3 + 46 x0x2x3 */

	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',0 }, {'z',1 }, {'z',3 } }, 15},
			{ { { 'x',1 }, {'y',2 }, {'z',3 } }, 78},
			{ { { 'z',1 }, {'x',3 },          }, -95},
			{ { { 'x',0 }, {'x',2 }, {'x',3 } }, 46},
		});

	Variable<default_qtraits> theta1;
	Variable<default_qtraits> theta2;
	Variable<default_qtraits> theta3;
	Variable<default_qtraits> theta4;
	Circuit<default_qtraits> cir;

	/*
	Circuit:
	  Rx(0, theta1) CNOT(0,1) Ry(0,theta2) CNOT(0,2) Rz(0,theta3)
	  Ry(1, theta2) --------- Rz(1,theta1) CNOT(1,3) Rx(1,theta4)
	  Ry(2, theta3) CNOT(2,3) Rz(2,theta1) --------- Ry(2,theta4)
	  Ry(3, theta1) --------- Rz(3,theta4) --------- Ry(3,theta2)	
	*/
	cir - RX(0, theta1)
		- RY(1, theta2)
		- RY(2, theta3)
		- RY(3, theta1)
		- CNOT(0, 1)
		- CNOT(2, 3)
		- RY(0, theta2)
		- RZ(1, theta1)
		- RZ(2, theta1)
		- RZ(3, theta4)
		- CNOT(0, 2)
		- CNOT(1, 3)
		- RZ(0, theta3)
		- RX(1, theta4)
		- RY(2, theta4)
		- RY(3, theta2);

	ObjectiveFunction obj({ theta1, theta2, theta3, theta4 }, cir, h);
	fp_t exp = -999;
	std::vector<fp_t> difs;
	std::vector<fp_t> compare_difs;

	/* theta1 = -0.123; theta2 = -0.899; theta3 = 1.111; theta4 = 2.001 */
	theta1.set_value(-0.123);
	theta2.set_value(-0.899);
	theta3.set_value(1.111);
	theta4.set_value(2.001);

	/* Find the expectation and differential at this point*/
	exp = obj.get_expectation_prob_noisefree();
	difs = obj.get_derivatives_prob_noisefree();

	fp_t delta = 1e-6; fp_t exp_plus = -999; fp_t exp_minus = -999;
	/* For every point, find the delta-diff */
	theta1.set_offset(delta);
	exp_plus = obj.get_expectation_prob_noisefree();
	theta1.set_offset(-2 * delta);
	exp_minus = obj.get_expectation_prob_noisefree();
	theta1.set_offset(delta);
	compare_difs.push_back((exp_plus - exp_minus) / 2 / delta);

	theta2.set_offset(delta);
	exp_plus = obj.get_expectation_prob_noisefree();
	theta2.set_offset(-2 * delta);
	exp_minus = obj.get_expectation_prob_noisefree();
	theta2.set_offset(delta);
	compare_difs.push_back((exp_plus - exp_minus) / 2 / delta);

	theta3.set_offset(delta);
	exp_plus = obj.get_expectation_prob_noisefree();
	theta3.set_offset(-2 * delta);
	exp_minus = obj.get_expectation_prob_noisefree();
	theta3.set_offset(delta);
	compare_difs.push_back((exp_plus - exp_minus) / 2 / delta);

	theta4.set_offset(delta);
	exp_plus = obj.get_expectation_prob_noisefree();
	theta4.set_offset(-2 * delta);
	exp_minus = obj.get_expectation_prob_noisefree();
	theta4.set_offset(delta);
	compare_difs.push_back((exp_plus - exp_minus) / 2 / delta);

	/* Expect Results:
		Exp =			1.46897
		Dif =			-22.5393, 4.20097, 47.8655, -53.0407
		Dif(compare) =  -22.5393, 4.20097, 47.8655, -53.0407	
	*/

	//std::cout << "Exp = " << exp << std::endl;
	//std::cout << "Dif = " << vec2str(difs, ", ") << std::endl;
	//std::cout << "Dif(compare) = " << vec2str(compare_difs, ", ") << std::endl;
	fp_t threshold = 1.e-4;

	res &= testcheck(roundn(1.46897, exp, threshold), "step 1");

	res &= testcheck(roundn(-22.5393, difs[0], threshold), "step 2");
	res &= testcheck(roundn(4.20097, difs[1], threshold), "step 3");
	res &= testcheck(roundn(47.8655, difs[2], threshold), "step 4");
	res &= testcheck(roundn(-53.0407, difs[3], threshold), "step 5");

	res &= testcheck(roundn(-22.5393, compare_difs[0], threshold), "step 6");
	res &= testcheck(roundn(4.20097, compare_difs[1], threshold), "step 7");
	res &= testcheck(roundn(47.8655, compare_difs[2], threshold), "step 8");
	res &= testcheck(roundn(-53.0407, compare_difs[3], threshold), "step 9");
	
	return res;
}

Test(DifferentialTest5) {
	TestResult res;

	using qid = typename default_qtraits::qidx_t;
	using fp_t = typename default_qtraits::value_t;

	/* H = 15 z0z1z3 + 78 x1y2z3 - 95 z1x3 + 46 x0x2x3 */

	Hamiltonian<default_qtraits> h(
		{
			{ { { 'z',0 }, {'z',1 }, {'z',3 } }, 15},
			{ { { 'x',1 }, {'y',2 }, {'z',3 } }, 78},
			{ { { 'z',1 }, {'x',3 },          }, -95},
			{ { { 'x',0 }, {'x',2 }, {'x',3 } }, 46},
		});

	Variable<default_qtraits> theta1;
	Variable<default_qtraits> theta2;
	Variable<default_qtraits> theta3;
	Variable<default_qtraits> theta4;
	Circuit<default_qtraits> cir;

	/*
	Circuit:
	  Rx(0, theta1) CNOT(0,1) Ry(0,theta2) CNOT(0,2) Rz(0,theta3)
	  Ry(1, theta2) --------- Rz(1,theta1) CNOT(1,3) Rx(1,theta4)
	  Ry(2, theta3) CNOT(2,3) Rz(2,theta1) --------- Ry(2,theta4)
	  Ry(3, theta1) --------- Rz(3,theta4) --------- Ry(3,theta2)
	*/
	cir - RX(0, theta1)
		- RY(1, theta2)
		- RY(2, theta3)
		- RY(3, theta1)
		- CNOT(0, 1)
		- CNOT(2, 3)
		- RY(0, theta2)
		- RZ(1, theta1)
		- RZ(2, theta1)
		- RZ(3, theta4)
		- CNOT(0, 2)
		- CNOT(1, 3)
		- RZ(0, theta3)
		- RX(1, theta4)
		- RY(2, theta4)
		- RY(3, theta2);

	ObjectiveFunction obj({ theta1, theta2, theta3, theta4 }, cir, h);
	fp_t exp = -999;
	std::vector<fp_t> difs;
	std::vector<fp_t> compare_difs;

	/* theta1 = -0.123; theta2 = -0.899; theta3 = 1.111; theta4 = 2.001 */
	theta1.set_value(-0.123);
	theta2.set_value(-0.899);
	theta3.set_value(1.111);
	theta4.set_value(2.001);

	size_t shots = 1000;
	/* Find the expectation and differential at this point*/

	/* Expect Results:
		Exp =			1.46897
		Dif =			-22.5393, 4.20097, 47.8655, -53.0407
	*/

	DefineTimer(timer);

	for (int i = 0; i < 1; ++i) {
		StartTimer(timer);
		exp = obj.get_expectation_montecarlo(shots, single_thread, noisy);
		difs = obj.get_derivatives_montecarlo(shots, single_thread, noisy);
		StopTimer(timer);

		std::cout << "Round " << i << ":" << std::endl;
		std::cout << "Exp = " << exp << std::endl;
		std::cout << "Dif = " << vec2str(difs, ", ") << std::endl;
	}

	for (int i = 0; i < 1; ++i) {
		StartTimer(timer);
		exp = obj.get_expectation_montecarlo(shots, multi_threads, noisy);
		difs = obj.get_derivatives_montecarlo(shots, multi_threads, noisy);
		StopTimer(timer);
		
		std::cout << "Round " << i << ":" << std::endl;
		std::cout << "Exp = " << exp << std::endl;
		std::cout << "Dif = " << vec2str(difs, ", ") << std::endl;
	}

	for (int i = 0; i < 1; ++i) {
		StartTimer(timer);
		exp = obj.get_expectation_montecarlo(shots, single_thread, noise_free);
		difs = obj.get_derivatives_montecarlo(shots, single_thread, noise_free);
		StopTimer(timer);

		std::cout << "Round " << i << ":" << std::endl;
		std::cout << "Exp = " << exp << std::endl;
		std::cout << "Dif = " << vec2str(difs, ", ") << std::endl;
	}

	for (int i = 0; i < 1; ++i) {
		StartTimer(timer);
		exp = obj.get_expectation_montecarlo(shots, multi_threads, noise_free);
		difs = obj.get_derivatives_montecarlo(shots, multi_threads, noise_free);
		StopTimer(timer);

		std::cout << "Round " << i << ":" << std::endl;
		std::cout << "Exp = " << exp << std::endl;
		std::cout << "Dif = " << vec2str(difs, ", ") << std::endl;
	}

	return res;
}