#pragma once
#include <global.h>
#include <testutil.h>
ns_easyquantum

/* Transform a decimal to binary (as string). The lowest digit to str[0]. */
std::string dec2bin(size_t n, size_t size);

/* transform a unsigned number to string. No safety gaurantee. */
template<size_t digit>
std::string to_binary_string(uint64_t n) {
	std::stringstream ss;
	ss << bitset<digit>(n);
	return ss.str();
}

/* Transform a binary sting to an integer */
size_t bin2dec(std::string s);

/* Compute 2^n */
template<typename uintv_t = uint64_t, typename uintd_t = uint32_t>
constexpr uintv_t pow2(uintd_t n) {
	return ((uintv_t)1) << (n);
}

template<typename uintv_t = uint64_t, typename uintd_t = uint32_t>
constexpr uintd_t log2(uintv_t n) {
	// get number of maximum output
	constexpr uintd_t ndigit = sizeof(n) * 8;
	if (n == 1) return 0;
	for (uintd_t i = 1; i < ndigit; ++i) {
		// 101 <= 1000 && 101 > 100 return 3
		if (n <= pow2<uintv_t, uintd_t>(i) && 
			n > pow2<uintv_t, uintd_t>(i - 1)) {
			return i - 1;
		}
	}
	return ndigit;
}

template<typename uintv_t = uint64_t>
constexpr uintv_t intsqrt(uintv_t n) {
	return (uintv_t)std::sqrt(n);
}

/* flip the binary number */
template<typename uintv_t = uint64_t, typename uintd_t = uint32_t>
constexpr uintv_t flip_bin(uintv_t n, uintd_t size) {
	return pow2<uintv_t>(size) - 1 - n;
}

/* reverse the binary number */
template<typename uintv_t = uint64_t, typename uintd_t = uint32_t>
uintv_t reverse_bin(uintv_t n, uintd_t size) {
	assert(n < ((uintv_t)1 << size), "Size not enough");

	std::string xs = dec2bin(n, size);
	std::reverse(xs.begin(), xs.end());

	return bin2dec(xs);
}

/* Assume real and imag are from a complex number*/
template<typename Ty>
Ty abs2(const Ty &real, const Ty &imag) {
	return real * real + imag * imag;
}

/* Complex multiplication without usage of complex type */
template<typename Ty>
constexpr void complex_multiplication(Ty& result_real, Ty& result_imag,
	const Ty real0, const Ty imag0, const Ty real1, const Ty imag1) {

	result_real = real0 * real1 - imag0 * imag1;
	result_imag = real0 * imag1 + real1 * imag0;
}

/* Return true if float number is e-close to 0. */
template<typename fp_t>
bool round0(fp_t val, fp_t e) {
	return abs(val) < abs(e);
}

/* Return true if two float number is e-close. */
template<typename fp_t>
bool roundn(fp_t val, fp_t n, fp_t e) {
	return round0(val - n, e);
}

/* math constant pi */
constexpr double pi = 3.1415926535897932384626433832795028841971693993751;
constexpr double PI = pi;
constexpr double Pi = pi;

/* math constant sqrt(2) */
constexpr double sq2 = 1.414213562373095048801688724209698078569671875377;
constexpr double SQ2 = sq2;

ns_end