#include <easy_quantum>

using_ns_easyquantum

/* Test Vectorization functions. */
Test(DefaultRandomEngineTest) {
	TestResult res;

	constexpr int seed = 100;
	DefaultRandomEngine rng1(seed);
	RandomEngine& rng1ref = rng1;
	std::vector<double> random_numbers11 = rng1ref(50);
	std::vector<double> random_numbers12 = rng1ref(50);

	DefaultRandomEngine rng2(seed + 1);
	RandomEngine& rng2ref = rng2;
	std::vector<double> random_numbers21 = rng2ref(50);
	std::vector<double> random_numbers22 = rng2ref(50);

	DefaultRandomEngine rng3(seed);
	RandomEngine& rng3ref = rng3;
	std::vector<double> random_numbers31 = rng3ref(50);
	std::vector<double> random_numbers32 = rng3ref(50);

	res &= testcheck(!is_same_container(random_numbers11, random_numbers12),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers21, random_numbers22),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers31, random_numbers32),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers21, random_numbers22), 
		"Different seeds generate same random numbers");
	res &= testcheck(is_same_container(random_numbers11, random_numbers31),
		"Same seeds generate different random numbers");
	res &= testcheck(is_same_container(random_numbers12, random_numbers32),
		"Same seeds generate different random numbers");

	return res;
}

/* Test Vectorization functions. */
Test(XCRandomEngineTest) {
	TestResult res;

	constexpr int seed = 100;
	XC_RandomEngine16807 rng1(seed);
	RandomEngine& rng1ref = rng1;
	std::vector<double> random_numbers11 = rng1ref(50);
	std::vector<double> random_numbers12 = rng1ref(50);

	XC_RandomEngine16807 rng2(seed + 1);
	RandomEngine& rng2ref = rng2;
	std::vector<double> random_numbers21 = rng2ref(50);
	std::vector<double> random_numbers22 = rng2ref(50);

	XC_RandomEngine16807 rng3(seed);
	RandomEngine& rng3ref = rng3;
	std::vector<double> random_numbers31 = rng3ref(50);
	std::vector<double> random_numbers32 = rng3ref(50);

	res &= testcheck(!is_same_container(random_numbers11, random_numbers12),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers21, random_numbers22),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers31, random_numbers32),
		"Generate same random numbers");
	res &= testcheck(!is_same_container(random_numbers21, random_numbers22),
		"Different seeds generate same random numbers");
	res &= testcheck(is_same_container(random_numbers11, random_numbers31),
		"Same seeds generate different random numbers");
	res &= testcheck(is_same_container(random_numbers12, random_numbers32),
		"Same seeds generate different random numbers");

	return res;
}