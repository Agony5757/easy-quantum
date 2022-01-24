#pragma once
#include <global.h>
#include <chrono>
#ifdef _WIN32
#include <Windows.h>
#endif
ns_easyquantum

struct red {
	static constexpr WORD attribute = FOREGROUND_RED | FOREGROUND_INTENSITY;
};
struct green {
	static constexpr WORD attribute = FOREGROUND_GREEN | FOREGROUND_INTENSITY;
};
struct blue {
	static constexpr WORD attribute = FOREGROUND_GREEN | FOREGROUND_INTENSITY;
};
struct yellow {
	static constexpr WORD attribute = FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY;
};
struct white {
	static constexpr WORD attribute = FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE | FOREGROUND_INTENSITY;
};

template<typename _Elem, typename _Traits, typename color_t, 
	std::enable_if_t<color_t::attribute == color_t::attribute, int> = 0
	>
std::basic_ostream<_Elem, _Traits>&
operator<<(std::basic_ostream<_Elem, _Traits>& s, color_t) {
#ifdef _WIN32
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hStdout, color_t::attribute);
	return s;
#else
	return s;
#endif
}

struct TestResult {
	std::optional<std::string> result;

	TestResult() {}
	TestResult(std::optional<std::string> res) :result(res) {}
	TestResult(const TestResult& res) { result = res.result; }
	TestResult& operator=(const TestResult& res) { result = res.result; }

	inline operator bool() {
		return !result.has_value();
	}

	inline operator std::optional<std::string>() {
		return result;
	}

	inline std::string to_string() {
		if (result.has_value()) return *result;
		return "";
	}

	inline TestResult& operator &= (TestResult& res) {
		if (!res.result.has_value())
			return *this;

		if (result.has_value()) {
			*result += " | ";
			*result += *res.result;
		}
		if (!result.has_value()) {
			result = res.result;
		}
		return *this;
	}

	inline TestResult operator &(TestResult& res) {
		TestResult newres(*this);
		newres &= res;
		return newres;
	}
};

using FuncType = std::function<TestResult()>;
using TestType = std::pair<FuncType, std::string>;

class TestSet {
private:
	TestSet() {}
	std::vector<TestType> tests;
public:
	static TestSet& getInstance() {
		static TestSet tests1;
		return tests1;
	}

	inline void run_all_tests() {
		std::cout << "***** Tests For easy_quantum *****" << std::endl;
		std::cout << "----------------------------------" << std::endl;
		std::cout << std::endl;
		int passed = 0;
		int tested = 0;
		std::vector<std::string> failed_tests;
		for (TestType &test : tests) {
			tested++;
			TestResult res = test.first();
			if (res.result.has_value()) {
				std::cout << red() << "FAIL " << white() << test.second << " | INFO: " << res.result.value() << std::endl;
				failed_tests.push_back(test.second);
			}
			else {
				passed++;
				std::cout << green() << "PASS " << white() << test.second << std::endl;
			}
		}
		std::cout << "Tests passed: " << passed << "/" << tested << "." << std::endl;
		if (passed == tested) {
			std::cout << "Grats! All tests passed." << std::endl;
		}
		else {
			std::cout << "The following tests failed: " << std::endl;
			for (auto t : failed_tests) {
				std::cout << t << std::endl;
			}
		}
	}

	void append_test(TestType &type) { tests.push_back(type); }
};

struct TestRegisterHelper {
	TestRegisterHelper(TestType &test) {
		TestSet::getInstance().append_test(test);
	}
};

#ifdef NO_ASSERT
constexpr bool no_assert = true;
#else
constexpr bool no_assert = false;
#endif

inline void assert(bool check_value, std::string info) {
	if (no_assert || check_value) {
		// pass
	}
	else {
		throw std::runtime_error(info.c_str());
	}
}

inline void assert(bool check_value) {
	if (no_assert || check_value) {
		// pass
	}
	else {
		throw std::runtime_error("Internal Error.");
	}
}

inline TestResult testcheck(bool check_value, std::string info) {
	if (!check_value) {
		return info;
	}
	else {
		return {};
	}
}

inline TestResult testcheck(bool check_value) {
	if (!check_value) {
		return "";
	}
	else {
		return {};
	}
}

ns_end

#define RegisterTest(test) static easyquantum::TestRegisterHelper _##test##_helper(TestType{&test,#test})

#define Test(testname) \
TestResult testname(); RegisterTest(testname); TestResult testname()

/* Time Elapsed */
#define DefineTimer(timer) std::chrono::steady_clock::time_point timer
#define StartTimer(timer) timer = std::chrono::steady_clock::now()
#define StopTimer_ms(timer) std::cout << "Elapsed: " \
							  << std::chrono::duration_cast<std::chrono::milliseconds> \
								(std::chrono::steady_clock::now() - timer).count()\
							  << " ms" << std::endl
#define StopTimer_s(timer) std::cout << "Elapsed: " \
							  << std::chrono::duration_cast<std::chrono::seconds> \
								(std::chrono::steady_clock::now() - timer).count()\
							  << " s" << std::endl
#define StopTimer(timer) StopTimer_ms(timer)

#define GetTimer_ms(timer) std::chrono::duration_cast<std::chrono::milliseconds> \
								(std::chrono::steady_clock::now() - timer).count()
#define GetTimer_s(timer) std::chrono::duration_cast<std::chrono::seconds> \
								(std::chrono::steady_clock::now() - timer).count()
#define GetTimer(timer) GetTimer_ms(timer)