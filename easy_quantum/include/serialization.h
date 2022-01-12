#pragma once
#include <global.h>

ns_easyquantum

/* Status: UNSTABLE. */

/* What's this:
	This class is for serialization of class/struct data to
	binary, then storing them to the disk. Also loading from
	disk and remaking these objects.
   How to use:
	Include this header. 1. From In your class, write two
	functions:
		 void to_binary(Binary &bin) const;
		 void from_binary(Binary &bin);
	Any integral type is already supported. If you have
	several	supported types of data: d1, d2, ..., you just
	compose	them by:
		bin & d1 & d2;

	And construct them by:
		bin >> d1;
		bin >> d2;

	Always remember the sequence.
*/

class Binary;

/* Status: UNSTABLE. */

/* What's this:
	This class is for serialization of class/struct data to
	binary, then storing them to the disk. Also loading from
	disk and remaking these objects.
   How to use:
	Include this header. 1. From In your class, write two
	functions:
		 void to_binary(Binary &bin) const;
		 void from_binary(Binary &bin);
	Any integral type is already supported. If you have
	several	supported types of data: d1, d2, ..., you just
	compose	them by:
		bin & d1 & d2;

	And construct them by:
		bin >> d1;
		bin >> d2;

	Always remember the sequence.
*/

class Binary;

template<typename Ty>
struct insert_impl {
	static Binary& insert(Binary& bin, const Ty& value);
};

template<typename Ty>
struct omit_impl {
	static Binary& omit(Binary& bin, Ty& value);
};

template<typename Ty>
struct insert_impl<std::vector<Ty>> {
	static Binary& insert(Binary& bin, const std::vector<Ty>& value);
};

template<typename Ty>
struct omit_impl<std::vector<Ty>> {
	static Binary& omit(Binary& bin, std::vector<Ty>& value);
};

template<typename Ty>
struct insert_impl<std::basic_string<Ty>> {
	static Binary& insert(Binary& bin, const std::basic_string<Ty>& value);
};

template<typename Ty>
struct omit_impl<std::basic_string<Ty>> {
	static Binary& omit(Binary& bin, std::basic_string<Ty>& value);
};

class Binary {
public:
	virtual void _insert(void* data, int size) = 0;
	virtual void* _omit(int size) = 0;
	virtual void _move_pointer(int size) = 0;
	virtual void omit_ready() = 0;

	template<typename Ty>
	Binary& operator&(const Ty& value) {
		return insert_impl<Ty>::insert(*this, value);
	};

	template<typename Ty>
	Binary& operator>>(Ty& value) {
		return omit_impl<Ty>::omit(*this, value);
	}
};

template<typename Ty>
Binary& insert_impl<Ty>::insert(Binary& bin, const Ty& value) {
	if constexpr (is_integral_v<Ty>) {
		int size = sizeof(int);
		bin._insert(reinterpret_cast<void*>(const_cast<Ty*>(&value)), size);
		return bin;
	}
	else {
		value.to_binary(bin);
	}
	return bin;
}

template<typename Ty>
Binary& omit_impl<Ty>::omit(Binary& bin, Ty& value) {
	if constexpr (is_integral_v<Ty>) {
		int size = sizeof(Ty);
		void* data = bin._omit(size);
		bin._move_pointer(size);
		value = *(reinterpret_cast<Ty*>(data));
	}
	else {
		value.from_binary(bin);
	}
	return bin;
}

template<typename Ty>
Binary& insert_impl<std::vector<Ty>>::insert(Binary& bin, const std::vector<Ty>& value) {
	bin.insert(value.size());
	for (decltype(value.size()) i = 0; i < value.size(); ++i) {
		bin.insert<Ty>(value[i]);
	}
	return bin;
}

template<typename Ty>
Binary &omit_impl<std::vector<Ty>>::omit(Binary &bin, std::vector<Ty>& value) {
	using size_type = decltype(value.size());
	size_type size;
	omit<size_type>(size);
	value.resize(size);
	for (size_type i = 0; i < size; ++i) {
		omit<Ty>(value[i]);
	}
	return bin;
}

template<typename Ty>
Binary& insert_impl<std::basic_string<Ty>>::insert(Binary &bin, const std::basic_string<Ty>& value) {
	auto size = value.size();
	auto pt = value.data();
	bin.insert(size);
	for (decltype(size) i = 0; i < size; ++i) {
		bin.insert<Ty>(pt[i]);
	}
	return bin;
}

template<typename Ty>
Binary &omit_impl<std::basic_string<Ty>>::omit(Binary &bin, std::basic_string<Ty>& value) {
	using size_type = decltype(value.size());
	size_type size;
	omit<size_type>(size);
	value.resize(size);
	for (size_type i = 0; i < size; ++i) {
		omit<Ty>(value[i]);
	}
	return bin;
}


template<typename chartype = char>
class BinaryVector : public Binary {
public:
	std::vector<chartype> x;
	typename std::vector<chartype>::size_type omit_iterator;

	BinaryVector() {
		omit_iterator = 0;
	}

	void omit_ready() {
		omit_iterator = 0;
	}

	void _move_pointer(int size) {
		omit_iterator += size;
	}

	void _insert(void* data, int size) {
		x.reserve(x.size() + size);
		const chartype* ptr = reinterpret_cast<const chartype*>(data);
		for (int i = 0; i < size; ++i) {
			x.push_back(ptr[i]);
		}
	}

	void* _omit(int size) {
		return x.data() + omit_iterator;
	}

	void write_to_file(std::basic_string<chartype> filename) {
		std::basic_ofstream<chartype> out(filename, std::ios::out);
		out << x.size();
		for (auto x_ : x) out << x_;
	}

	void write_to_binary_file(std::basic_string<chartype> filename) {
		std::basic_ofstream<chartype> out(filename, std::ios::out | std::ios::binary);
		typename decltype(x)::size_type size = x.size();
		out.write(reinterpret_cast<const chartype*>(&size), sizeof(size));
		for (auto x_ : x) out.put(x_);
	}

	void load_from_binary_file(std::basic_string<chartype> filename) {
		std::basic_ifstream<chartype> in(filename, std::ios::in | std::ios::binary);
		x.clear();
		chartype ptsize[128];
		using size_type = typename decltype(x)::size_type;
		in.read(ptsize, sizeof(size_type));
		size_type size = *reinterpret_cast<size_type*>(ptsize);
		for (size_type i = 0; i < size; ++i) {
			x.push_back(in.get());
		}
	}

	void test_output() {
		std::cout << "BIN: ";
		for (auto i : x) {
			std::cout << i << " ";
		}
		std::cout << std::endl;
	}
};

ns_end

