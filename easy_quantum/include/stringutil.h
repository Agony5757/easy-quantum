#pragma once
#include <global.h>

ns_easyquantum

/* Merge two strings */
template<typename char_t>
std::basic_string<char_t> merge_strings(
	std::vector<std::basic_string<char_t>> strings, 
	std::basic_string<char_t> delim) {

	std::basic_string<char_t> final_string;
	for (decltype(strings.size()) i = 0;
		i < strings.size() - 1; ++i) {
		final_string += strings[i];
		final_string += delim;
	}
	final_string += strings.back();
	return final_string;
}

/* String split */
template<typename char_t>
std::vector<std::basic_string<char_t>>
split(const std::basic_string<char_t>& s, const std::basic_string<char_t>& delimiters) {

	std::vector<std::basic_string<char_t>> tokens;
	std::basic_string<char_t>::size_type lastPos = s.find_first_not_of(delimiters, 0);
	std::basic_string<char_t>::size_type pos = s.find_first_of(delimiters, lastPos);
	while (std::basic_string<char_t>::npos != pos || 
		std::basic_string<char_t>::npos != lastPos) {

		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delimiters, pos);
		pos = s.find_first_of(delimiters, lastPos);
	}
	return tokens;
}

/* Trim: remove all blank characters at the head/tail. */
std::wstring trim(std::wstring str, const std::wstring &blank_chars = L" ");

/* Trim: remove all blank characters at the head/tail. */
std::string trim(std::string str, const std::string &blank_chars = " ");

/* Use regular expression to check if a string is integer */
template<typename char_t>
bool is_integer_str(const std::basic_string<char_t> &str) {
	std::basic_regex<char_t> re("[0-9]+");
	return std::regex_match(str, re);
}

template<>
inline bool is_integer_str<wchar_t>(const std::basic_string<wchar_t> &str) {
	std::basic_regex<wchar_t> re(L"[0-9]+");
	return std::regex_match(str, re);
}

#define SHOWPOINT true
#define NOSHOWPOINT false

/* Transform a vector to a string */
template<typename Ty>
std::string vec2str(std::vector<Ty> vec, std::string sep = ", ", bool showpt = NOSHOWPOINT) {
	if (vec.size() == 0) {
		return std::string();
	}
	std::stringstream ss;
	if (showpt)
		ss << std::showpoint;

	for (auto iter = vec.begin(); iter != (vec.end() - 1); ++iter) {
		ss << *iter << sep;
	}
	ss << vec.back();
	return ss.str();
}

/* Transform a vector to a string */
template<typename ContainerTy, typename Ty = ContainerTy::value_type>
std::string vec2str(ContainerTy vec, std::function<std::string(const Ty&)> func,
	std::string sep = ", ", bool showpt = NOSHOWPOINT) {

	std::stringstream ss;
	if (showpt)
		ss << std::showpoint;

	if (vec.size() == 0) {
		return "";
	}

	for (auto iter = vec.begin(); iter != std::prev(vec.end()); ++iter) {
		ss << func(*iter) << sep;
	}
	ss << func(*std::prev(vec.end()));
	return ss.str();
}

/* Transform a vector to a wstring */
template<typename Ty>
std::wstring vec2str(std::vector<Ty> vec, std::wstring sep = L", ", bool showpt = NOSHOWPOINT) {
	std::wstringstream ss;
	if (showpt)
		ss << std::showpoint;

	for (auto iter = vec.begin(); iter != (vec.end() - 1); ++iter) {
		ss << *iter << sep;
	}
	ss << vec.back();
	return ss.str();
}

template<typename Ty>
std::string complex_to_string(const std::complex<Ty>& c) {
	std::stringstream ss;
	ss << std::to_string(c.real);
	if (c.imag >= 0) {
		ss << " + ";
	}
	else {
		ss << " ";
	}
	ss << c.imag << "i";
	return ss.str();
}

ns_end