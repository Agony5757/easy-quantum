#include <easy_quantum>

using_ns_easyquantum

/* Test bin2dec and dec2bin*/
Test(MathUtilTest1) {
	TestResult res;

	size_t x = 1022;
	std::string xs1 = dec2bin(x, 10);
	std::string xs2 = dec2bin(x, 12);

	size_t x1 = bin2dec(xs1);
	size_t x2 = bin2dec(xs2);

	res &= testcheck(x1 == x2, "Bad Result.");

	return res;
}

/* Test reverse_bin */
Test(MathUtilTest2) {
	TestResult res;
	size_t x = 1022;
	unsigned int size = 12;
	std::string xs = dec2bin(x, size);
	
	reverse(xs.begin(), xs.end());
	size_t x1 = bin2dec(xs);

	size_t xrev = reverse_bin(x, size);
	std::string xrevs = dec2bin(xrev, size);
	/*char buf[1024];
	sprintf_s(buf, "Expect: %ld == %ld | Binary: %s, %s", 
		x1, xrev, xs.c_str(), xrevs.c_str());*/
	res &= testcheck(xrev == x1);

	return res;
}

/* Test flip_bin */
Test(MathUtilTest3) {
	TestResult res;
	size_t x = 1022;
	unsigned int size = 12;
	std::string xs = dec2bin(x, size);

	for_all(xs, [](char& c) {
		if (c == '0') c = '1';
		else if (c == '1') c = '0';
		else assert(false);
	});
	size_t x1 = bin2dec(xs);

	size_t xrev = flip_bin(x, size);
	std::string xrevs = dec2bin(xrev, size);
	/*char buf[1024];
	sprintf_s(buf, "Expect: %ld == %ld | Binary: %s, %s",
		x1, xrev, xs.c_str(), xrevs.c_str());*/
	res &= testcheck(xrev == x1);

	return res;
}
