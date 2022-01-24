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
