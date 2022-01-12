#include <mathutil.h>

ns_easyquantum

std::string dec2bin(size_t n, size_t size) {
	assert(n < ((size_t)1 << size), "Size not enough");
	std::string binstr = "";
	for (size_t i = 0; i < size; ++i) {
		binstr = (char)((n & 1) + '0') + binstr;
		n >>= 1;
	}
	return binstr;
}

size_t bin2dec(std::string s) {
	assert(s.size() <= 64, "Too long string.");
	size_t ret = 0;
	for (size_t i = 0; i < s.size(); ++i) {
		size_t real_i = s.size() - i - 1;
		if (s[i] == '1') {
			ret += (1ull << real_i);
		}
		else if (s[i] == '0') {
			ret += 0;
		}
		else assert(false, "Bad binary string");
	}
	return ret;
}

ns_end