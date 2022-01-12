#include <stringutil.h>

ns_easyquantum

std::string merge_strings(std::vector<std::string> strings) {
	std::string final_string;
	for (auto _string : strings) {
		final_string += _string;
		final_string += " ";
	}
	return final_string;
}

std::vector<std::string> split(const std::string& str, const std::string& delim) {
	std::vector<std::string> res;
	if ("" == str) return res;
	char *strs = new char[str.length() + 1];
	strcpy(strs, str.c_str());

	char *d = new char[delim.length() + 1];
	strcpy(d, delim.c_str());

	char *p = strtok(strs, d);
	while (p) {
		std::string s = p;
		res.push_back(s);
		p = strtok(NULL, d);
	}
	delete[] d;
	delete[] strs;

	return res;
}

std::string trim(std::string s, const std::string &blank_chars) {
	if (s.empty())
		return s;
	s.erase(0, s.find_first_not_of(blank_chars));
	s.erase(s.find_last_not_of(blank_chars) + 1);
	return s;
}

std::wstring trim(std::wstring s, const std::wstring &blank_chars) {
	if (s.empty())
		return s;
	s.erase(0, s.find_first_not_of(blank_chars));
	s.erase(s.find_last_not_of(blank_chars) + 1);
	return s;
}

ns_end