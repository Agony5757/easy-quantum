#include <vectorutil.h>

using_ns_easyquantum

/* Test Vectorization functions. */
Test(VectorUtilTest1) {
	TestResult res;

	std::vector<int> a = { 1,2,3 };
	std::vector<int> b = { 2,4,6 };	
	int c = 2;

	auto d = a + b;
	auto e = a * c;

	res &= testcheck(d[0] == 3);
	res &= testcheck(d[1] == 6);
	res &= testcheck(d[2] == 9);
	res &= testcheck(e[0] == 2);
	res &= testcheck(e[1] == 4);
	res &= testcheck(e[2] == 6);

	return res;
}

/* Test is_in */
Test(VectorUtilTest2) {
	struct A {
		int a;
		inline bool operator==(const A& a_) const {
			return a == a_.a;
		}
		A(int a_) :a(a_) {}
	};
	std::vector<A> as = { 1,2,3,4,5 };

	TestResult res;
	res &= testcheck(is_in(as, A(1)), "Step1/1");
	res &= testcheck(is_in(as, A(5)), "Step1/2");
	res &= testcheck(!is_in(as, A(6)), "Step1/3");
	res &= testcheck(!is_in(as, A(-1)), "Step1/4");

	res &= testcheck(is_in(as, 2, 
		[](const A& a, const int b) {return a.a == b; }), "Step2/1");

	res &= testcheck(is_in(as, 3,
		[](const A& a, const int b) {return a.a == b; }), "Step2/2");

	res &= testcheck(!is_in(as, 6,
		[](const A& a, const int b) {return a.a == b; }), "Step2/3");

	return res;
}

/* Test erase & is_same_container */
Test(VectorUtilTest3) {
	TestResult res;

	std::vector<int> a = { 1,2,3,4,5,6,7,8,9,10 };
	std::vector<int> b = { 1,2,9,10 };
	std::vector<int> c_ = { 3,4,5,6,7,8 };
	std::vector<int> c = erase_many(a, b);

	res &= testcheck(is_same_container(c_, c), "Step1/1");
	res &= testcheck(!is_same_container(a, c), "Step1/2");
	res &= testcheck(!is_same_container(b, c), "Step1/3");

	return res;
}

/* Test for_all */
Test(VectorUtilTest4) {
	TestResult res;

	std::vector<int> a = { 1,2,3,4,5 };
	std::vector<int> a_ = a;
	for_all(a, [](int& i) { i++; });
	std::vector<int> b = { 2,3,4,5,6 };
	res &= testcheck(is_same_container(a, b), "Step1/1");
	res &= testcheck(is_same_container(a, a), "Step1/2");
	res &= testcheck(!is_same_container(a_, a), "Step1/3");

	return res;
}

/* Test for_all */
Test(VectorUtilTest5) {
	TestResult res;

	std::vector<int> a = { 1,2,3,4,5 };
	std::vector<int> a_ = a;
	std::vector<int> b = { 5,4,3,2,1 };
	std::vector<int> c = { 1,2,3,4,5,5,4,3,2,1 };
	merge_into(a, b);

	res &= testcheck(is_same_container(a, c), "Step1/1");
	res &= testcheck(!is_same_container(a, b), "Step1/2");
	res &= testcheck(!is_same_container(a_, a), "Step1/3");

	return res;
}