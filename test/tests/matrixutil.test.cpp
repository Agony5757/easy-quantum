#include <easy_quantum>

using_ns_easyquantum

Test(MatrixUtilTest1) {
	TestResult res;

	using matrix_t = DenseMatrix<double>;
	using vector_t = DenseVector<double>;

	matrix_t M1 = randmat<double>(10);
	vector_t V1 = randvec<double>(10);

	matrix_t M2 = randmat<double>(10);
	vector_t V2 = randvec<double>(10);

	res &= testcheck(M1 * V1 == M1 * V1, "Step 1/2");
	res &= testcheck(M1 * V1 != M2 * V2, "Step 2/2");

	return res;
}

Test(MatrixUtilTest2) {
	TestResult res;

	using matrix_t = DenseMatrix<double>;
	using vector_t = DenseVector<double>;

	vector_t v = { 1,2,3,4,5 };
	res &= testcheck(v.size == 5, "Step 1");
	res &= testcheck(
		v[0] == 1 &&
		v[1] == 2 &&
		v[2] == 3 &&
		v[3] == 4 &&
		v[4] == 5, "Step 2");

	vector_t v0 = { 5,4,3,2,1 };

	vector_t v1 = v0 + v;

	res &= testcheck(
		v1[0] == 6 &&
		v1[1] == 6 &&
		v1[2] == 6 &&
		v1[3] == 6 &&
		v1[4] == 6, "Step 3");

	return res;
}

Test(MatrixUtilTest3) {
	TestResult res;

	using matrix_t = DenseMatrix<double>;
	using vector_t = DenseVector<double>;

	vector_t v{ 1,-1 };
	matrix_t m{ 
		2,3,
		4,5 
	};

	matrix_t m0{
		-1,1,
		1,-1
	};

	vector_t v1{ -1,-1 };
	matrix_t m1{
		1, -1,
		1, -1
	};

	assert(v1 == m * v, "Step 1/2");
	assert(m1 == m * m0, "Step 2/2");

	return res;
}