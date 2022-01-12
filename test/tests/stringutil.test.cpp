#include <easy_quantum>

using_ns_easyquantum

// Test merge string
Test(StringUtilTest1) {
	TestResult res;
	
	std::vector<std::string> strings = {
		"str1", "str2", "str3", "str4"
	};
	
	std::string str_= "str1, str2, str3, str4";
	std::string str = merge_strings(strings, std::string(", "));
	
	res &= testcheck(str_ == str);

	return res;
}

// Test split string
Test(StringUtilTest2) {
	TestResult res;

	{
		std::string str = "str1,,str2, ,str3,  ,str4,,,";
		std::vector<std::string> strs = split(str, std::string(","));
		
		std::vector<std::string> strs_ = {
			"str1","str2"," ","str3","  ","str4"
		};
		res &= testcheck(is_same_container(strs_, strs), "Step 1/1");
	}
	{
		std::string str = "str1,,str2, ,str3,  ,str4,,,";
		std::vector<std::string> strs = split(str, std::string(", "));

		std::vector<std::string> strs_ = {
			"str1","str2","str3","str4"
		};
		res &= testcheck(is_same_container(strs_, strs), "Step 1/2");
	}

	return res;
}

// Test trim
Test(StringUtilTest3) {
	TestResult res;
	std::string str = "  AAA ,  AAA   ,  ";
	std::string a = trim(str, ", ");
	std::string b = trim(str, ",");
	std::string c = trim(str, " ");
	
	res &= testcheck(a == "AAA ,  AAA");
	res &= testcheck(b == str);
	res &= testcheck(c == "AAA ,  AAA   ,");

	return res;
}

// Test trim
Test(StringUtilTest4) {
	TestResult res;
	std::string str = "  AAA ,  AAA   ,  ";
	std::string a = trim(str, ", ");
	std::string b = trim(str, ",");
	std::string c = trim(str, " ");

	res &= testcheck(a == "AAA ,  AAA");
	res &= testcheck(b == str);
	res &= testcheck(c == "AAA ,  AAA   ,");

	return res;
}