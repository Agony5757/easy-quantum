#pragma once

#include <global.h>
#include <memoryutil.h>
#include <randomutil.h>

ns_easyquantum

template<typename Ty, typename uint, typename allocator> 
struct DenseVector;

/* Class of dense matrix. */
template<typename Ty, 
	typename uint = unsigned int, 
	typename allocator = safe_allocator<Ty>>
struct DenseMatrix {

	using size_type = uint;
	using uint_t = uint;

	size_type size = 0;
	Ty* data = nullptr;

	DenseMatrix() {}

	void initialize(size_type size_, allocator A = allocator()) {
		assert(size == 0 && data == nullptr, "Only empty matrix can be initialized.");

		size = size_;
		// data = new Ty[size * size];
		data = A.allocate(size * size);
		
		for (size_type i = 0; i < size * size; ++i) data[i] = 0;
	}

	DenseMatrix(std::initializer_list<Ty> list) {
		size_type s = (size_type)list.size();
		size_type size0 = intsqrt(s);
		assert(size0 * size0 == s, "Input size should be n*n.");
		initialize(size0);
		auto p = list.begin();
		for (size_type i = 0; i < size0 * size0; ++i) {
			data[i] = *p;
			++p;
		}
	}

	/* Create a size*size matrix, all values cleared to 0. */
	DenseMatrix(size_type size_, allocator A = allocator()) {
		initialize(size_, A);
	}

	/* Copy constructor. Deep copy. */
	DenseMatrix(const DenseMatrix<Ty, uint_t, allocator>& m) {
		size = m.size;
		// data = new Ty[m.size*m.size];
		allocator A;
		data = A.allocate(m.size * m.size);

		for (size_type i = 0; i < size*size; ++i) data[i] = m.data[i];
	}

	/* Get a reference to one element. You can modify this element. */
	Ty& get(size_type x, size_type y) {
		return data[x * size + y];
	}

	/* Same as get(x,y). */
	Ty& operator()(size_type x, size_type y) {
		return get(x, y);
	}

	const Ty& at(size_type x, size_type y) const {
		return data[x * size + y];
	}

	/* Matrix multiplication without optimization. */
	DenseMatrix<Ty, uint_t, allocator> operator*(const DenseMatrix<Ty, uint_t, allocator>& m) const {
		DenseMatrix<Ty, uint_t, allocator> newm(size);
		for (size_type i = 0; i < size; ++i) {
			for (size_type j = 0; j < size; ++j) {
				for (size_type k = 0; k < size; ++k) {
					newm(i, j) += (*this).at(i, k) * m.at(k, j);
				}
			}
		}
		return newm;
	}

	/* Stringify this matrix. */
	std::string to_string() const {
		std::stringstream ss;
		for (size_type i = 0; i < size; ++i) {
			for (size_type j = 0; j < size; ++j) {
				ss << at(i, j) << " ";
			}
			ss << endl;
		}
		return ss.str();
	}

	/* Write the matrix into a matlab readable file. */
	void write_matlab_file(std::string filename) const {
		ofstream out(filename, ios::out);
		for (size_type i = 0; i < size; ++i) {
			for (size_type j = 0; j < size; ++j) {
				out << at(i, j) << " ";
			}
			out << endl;
		}
	}

	/* Destructor. Delete all allocated spaces. */
	~DenseMatrix() {
		delete[] data;
	}
};

/* Class of dense vector. */
template<typename Ty, typename uint = unsigned int, typename allocator = safe_allocator<Ty>>
struct DenseVector {
	using size_type = uint;

	size_type size = 0;
	Ty* data = nullptr;

	DenseVector() {}

	void initialize(size_type size_, allocator A = allocator()) {
		assert(size == 0 && data == nullptr, "Only empty matrix can be initialized.");
		size = size_;
		// data = new Ty[size];
		data = A.allocate(size);
		for (size_type i = 0; i < size; ++i) data[i] = 0;
	}

	void clear() {
		for (size_type i = 0; i < size; ++i) data[i] = 0;
	}

	/* Create a size*size matrix, all values cleared to 0. */
	DenseVector(size_type size_, allocator A = allocator()) {
		initialize(size_, A);
	}

	DenseVector(std::initializer_list<Ty> list) {
		initialize((size_type)list.size());
		size_type i = 0;
		for (Ty it : list) {
			data[i] = it;
			++i;
		}
	}

	/* Copy Constructor. Deep copy. */
	DenseVector(const DenseVector<Ty, uint, allocator>& m) {
		size = m.size;
		// data = new Ty[m.size];
		allocator A;
		data = A.allocate(m.size);
		for (size_type i = 0; i < size; ++i) data[i] = m.data[i];
	}
	
	/* Assign operator. First free the target space, then deep copy. */
	DenseVector<Ty, uint, allocator>& operator=(const DenseVector<Ty, uint, allocator>& m) {
		allocator A;
		A.deallocate(data, -1);
		size = m.size;
		data = A.allocate(m.size);
		for (size_type i = 0; i < size; ++i) data[i] = m.data[i];
		return *this;
	}

	/* Vector addition operator. */
	DenseVector<Ty, uint, allocator> operator+(const DenseVector<Ty, uint, allocator>& v) const {
		assert(v.size == size);
		DenseVector<Ty, uint, allocator> vout(*this);
		for (size_type i = 0; i < size; ++i)
			vout.data[i] += v.data[i];
		return vout;
	}

	/* Vector subtraction operator. */
	DenseVector<Ty, uint, allocator> operator-(const DenseVector<Ty, uint, allocator>& v) const {
		assert(v.size == size);
		DenseVector<Ty, uint, allocator> vout(*this);
		for (size_type i = 0; i < size; ++i)
			vout.data[i] -= v.data[i];
		return vout;
	}

	/* Get a reference to one element. You can modify this element. */
	Ty& get(size_type x) {
		assert(x >= 0 && x < size);
		return data[x];
	}

	/* Same as get(x). */
	Ty& operator()(size_type x) {
		return get(x);
	}

	/* Same as get(x). */
	Ty& operator[](size_type x) {
		return get(x);
	}

	const Ty& at(size_type x) const {
		return data[x];
	}

	/* Compute the l2 norm. */
	Ty norm2() {
		Ty sum = 0;
		for (size_type i = 0; i < size; ++i) sum += (data[i] * data[i]);
		return sqrt(sum);
	}

	/* Compute the l-inf norm. */
	Ty normInf() {
		Ty ninf = 0;
		for (size_type i = 0; i < size; ++i) {
			if (abs(data[i]) > ninf) {
				ninf = abs(data[i]);
			}
		}
		return ninf;
	}

	/* Stringify this vector. */
	std::string to_string() {
		std::stringstream ss;
		for (size_type i = 0; i < size; ++i) {
			ss << data[i] << endl;
		}
		return ss.str();
	}

	/* Write the vector into a matlab readable file. */
	void write_matlab_file(std::string filename) {
		ofstream out(filename, ios::out);
		for (size_type i = 0; i < size - 1; ++i) {
			out << get(i) << ";";
		}
		out << get(size - 1);
	}

	/* Compute the maximum value. 
	Output via the argument, return the index. 
	It will only return the first maximum (if there are more than 1 maximum). 
	*/
	size_type compute_max(Ty& maxvalue) {
		size_type idx = 0;
		maxvalue = get(0);
		for (size_type i = 1; i < size; ++i) {
			if (get(i) > maxvalue) {
				maxvalue = get(i);
				idx = i;
			}
		}
		return idx;
	}

	/* Compute the maximum value of its absolute value.
	Output via the argument, return the index.
	It will only return the first maximum (if there are more than 1 maximum).
	*/
	size_type maxabs(Ty& maxvalue) {
		size_type idx = 0;
		maxvalue = abs(get(0));
		for (size_type i = 1; i < size; ++i) {
			if (abs(get(i)) > maxvalue) {
				maxvalue = abs(get(i));
				idx = i;
			}
		}
		return idx;
	}

	/* Destructor. Delete all allocated spaces. */
	~DenseVector() {
		delete[] data;
	}
};

/* Matrix multipies vector. Return a new vector. */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> operator*(const DenseMatrix<Ty, uint_t, allocator_t> &m, 
	const DenseVector<Ty, uint_t, allocator_t> &v) {
	assert(m.size == v.size, "Bad size");
	DenseVector<Ty, uint_t, allocator_t> v2(m.size);
	//#pragma omp parallel
	for (DenseVector<Ty, uint_t, allocator_t>::size_type i = 0; i < m.size; ++i) {
		for (DenseVector<Ty, uint_t, allocator_t>::size_type j = 0; j < m.size; ++j) {
			v2(i) += m.at(i, j) * v.at(j);
		}
	}
	return v2;
}

/* Matrix multipies constant. */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> operator*(
	const DenseVector<Ty, uint_t, allocator_t> &v, const Ty& a) {
	DenseVector<Ty, uint_t, allocator_t> v2(v.size);
	//#pragma omp parallel
	for (DenseVector<Ty, uint_t, allocator_t>::size_type i = 0; i < v.size; ++i) {
		v2(i) = v.at(i) * a;
	}
	return v2;
}

/* Matrix divides vector. */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> operator/(
	const DenseVector<Ty, uint_t, allocator_t> v, const Ty& a) {
	return v * Ty(1.0 / a);
}

/* Get a column of the matrix. Return a new vector (not reference). */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> get_column(const DenseMatrix<Ty, uint_t, allocator_t> &A, 
	const typename DenseMatrix<Ty, uint_t, allocator_t>::size_type col) {

	DenseVector<Ty, uint_t, allocator_t> v = DenseVector<Ty, uint_t, allocator_t>(A.size);
	for (DenseVector<Ty, uint_t, allocator_t>::size_type i = 0; i < A.size; ++i) {
		v(i) = A(i, col);
	}
	return v;
}

constexpr double default_tolerance = 1e-3;

template<typename Ty, typename uint_t, typename allocator_t>
bool compare_vec_l2(const DenseVector<Ty, uint_t, allocator_t> &v1, 
	const DenseVector<Ty, uint_t, allocator_t> &v2, Ty tolerance = default_tolerance) {

	assert(v1.size == v2.size, "Bad size");
	assert(tolerance > 0, "Bad input.");
	if ((v1 - v2).norm2() < tolerance) return true;
	return false;
}

/* Compare two vectors. Tolerance = 1e-3.
Return true if the difference of two vectors has norm2 less than the tolerance.*/
template<typename Ty, typename uint_t, typename allocator_t>
bool operator==(const DenseVector<Ty, uint_t, allocator_t> &v1, 
	const DenseVector<Ty, uint_t, allocator_t> &v2) {
	return compare_vec_l2(v1, v2, default_tolerance);
}

/* Compare two vectors. Tolerance = 1e-3.
Return false if the difference of two vectors has norm2 less than the tolerance.*/
template<typename Ty, typename uint_t, typename allocator_t>
bool operator!=(const DenseVector<Ty, uint_t, allocator_t> &v1, 
	const DenseVector<Ty, uint_t, allocator_t> &v2) {
	return !(v1 == v2);
}

template<typename Ty, typename uint_t, typename allocator_t>
bool compare_two_matrices(const DenseMatrix<Ty, uint_t, allocator_t> &m1,
	const DenseMatrix<Ty, uint_t, allocator_t> &m2, const Ty tolerance = default_tolerance) {
	assert(m1.size == m2.size, "Bad size.");
	assert(tolerance > 0, "Bad input.");
	for (DenseMatrix<Ty, uint_t, allocator_t>::size_type i = 0; i < m1.size * m1.size; ++i) {
		if (abs(m1.data[i] - m2.data[i]) > tolerance) return false;
	}
	return true;
}

/* Compare two matrices. Tolerance = 1e-3.*/
template<typename Ty, typename uint_t, typename allocator_t>
bool operator==(const DenseMatrix<Ty, uint_t, allocator_t> &m1, 
	const DenseMatrix<Ty, uint_t, allocator_t> &m2) {
	return compare_two_matrices(m1, m2);
}

/* Compare two matrices. Tolerance = 1e-3.*/
template<typename Ty, typename uint_t, typename allocator_t>
bool operator!=(const DenseMatrix<Ty, uint_t, allocator_t> &m1, 
	const DenseMatrix<Ty, uint_t, allocator_t> &m2) {
	return !(m1==m2);
}

/* Compare two vectors. Print the different elements. */
template<typename Ty, typename uint_t, typename allocator_t>
void check_vec(const DenseVector<Ty, uint_t, allocator_t> &v1, 
	const DenseVector<Ty, uint_t, allocator_t> &v2, const Ty tolerance = default_tolerance) {
	assert(v1.size == v2.size, "Bad size");
	assert(tolerance > 0, "Bad input.");
	for (DenseVector<Ty, uint_t, allocator_t>::size_type i = 0; i < v1.size; ++i) {
		if (abs(v1(i) - v2(i)) > tolerance) {
			cout << "Bad: " << i << "\t" << "v1: " << v1(i) << "\tv2:" << v2(i) << endl;
		}
	}
}

/* Swap two rows in both a matrix and a vector (used in Gaussian elimination). */
template<typename Ty, typename uint_t, typename allocator_t>
void swap_two_rows(DenseMatrix<Ty, uint_t, allocator_t>& A, DenseVector<Ty, uint_t, allocator_t>& b, 
	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type row1, typename DenseMatrix<Ty, uint_t, allocator_t>::size_type row2) {
	swap(b(row1), b(row2));
	//#pragma omp parallel
	for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type i = 0; i < A.size; ++i) {
		swap(A(row1, i), A(row2, i));
	}
}

/* Perform row elimination step (used in Gaussian elimination). */
template<typename Ty, typename uint_t, typename allocator_t>
void row_elimination(DenseMatrix<Ty, uint_t, allocator_t>& A, DenseVector<Ty, uint_t, allocator_t>& b, 
	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type row, typename DenseMatrix<Ty, uint_t, allocator_t>::size_type row2) {
	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type size = A.size;
	assert(A(row, row) != 0 && A.size == b.size, "Bad Gaussian Solver.");
	if (A(row2, row) == 0) return;
	Ty coef = A(row2, row) / A(row, row);
	for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j = row; j < size; ++j) {
		A(row2, j) -= (A(row, j)*coef);
	}
	b(row2) -= (b(row)*coef);
}

/* Perform column elimination step (used in Gaussian elimination). */
template<typename Ty, typename uint_t, typename allocator_t>
void column_elimination(DenseMatrix<Ty, uint_t, allocator_t>& A, DenseVector<Ty, uint_t, allocator_t> &b, DenseVector<Ty, uint_t, allocator_t> &x) {
	static_assert(is_same_v<typename DenseMatrix<Ty, uint_t, allocator_t>::size_type, DenseVector<Ty, uint_t, allocator_t>::size_type>, "Bad type.");

	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type size = A.size;
	x(size - 1) = b(size - 1) / A(size - 1, size - 1);
	for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type i = size - 2; i >= 0; --i) {
		Ty b0 = b(i);
		for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j = i + 1; j < size; ++j) {
			b0 -= (A(i, j)*x(j));
		}
		x(i) = b0 / A(i, i);
	}
}

/* Linear solver using gaussian elimination. Return the result vector. */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> gaussian_linear_solver(DenseMatrix<Ty, uint_t, allocator_t>& A, DenseVector<Ty, uint_t, allocator_t> &b) {
	static_assert(is_same_v<
		typename DenseMatrix<Ty, uint_t, allocator_t>::size_type, typename DenseVector<Ty, uint_t, allocator_t>::size_type>, 
		"Bad type.");
	
	assert(A.size == b.size, "Bad size");
	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type n = A.size;
	DenseVector<Ty, uint_t, allocator_t> x(n);

	for (int i = 0; i < n; ++i) {
		// find maximum row
		Ty maximum = 0.0;
		typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j_ = i;
		for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j = i + 1; j < n; ++j) {
			if (abs(A(j, i)) > maximum) {
				j_ = j;
				maximum = abs(A(j, i));
			}
		}
		if (A(j_, i) == 0) {
			DenseVector<Ty, uint_t, allocator_t> vi = get_column(A, i);
			cout << vi.norm2();
			assert(false, "No solution.");
		}
		for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j = i + 1; j < n; ++j) {
			row_elimination(A, b, i, j);
		}
	}
	column_elimination(A, b, x);
	return x;
}

/* Created a random matrix. */
template<typename Ty, typename uint_t = unsigned int, typename allocator_t = safe_allocator<Ty>>
DenseMatrix<Ty, uint_t, allocator_t> randmat(
	typename DenseMatrix<Ty, uint_t, allocator_t>::size_type size, 
	RandomEngine* rng = nullptr) {
	RandomEngine *ptr = rng;
	
	bool allocated = false;
	if (rng == nullptr) {
		ptr = new DefaultRandomEngine();
		allocated = true;
	}

	RandomEngine& rand = *ptr;

	DenseMatrix<Ty, uint_t, allocator_t> m(size);
	for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type i = 0; i < size; ++i) {
		for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type j = 0; j < size; ++j) {
			m(i, j) = rand();
		}
	}
	if (allocated) {
		delete ptr;
	}

	return m;
}

/* Create a random vector. */
template<typename Ty, typename uint_t = unsigned int, typename allocator_t = safe_allocator<Ty>>
DenseVector<Ty, uint_t, allocator_t> randvec(
	typename DenseVector<Ty, uint_t, allocator_t>::size_type size, 
	RandomEngine *rng = nullptr) {

	RandomEngine* ptr = rng;

	bool allocated = false;
	if (rng == nullptr) {
		ptr = new DefaultRandomEngine();
		allocated = true;
	}

	RandomEngine& rand = *ptr;

	DenseVector<Ty, uint_t, allocator_t> v(size);
	for (typename DenseMatrix<Ty, uint_t, allocator_t>::size_type i = 0; i < size; ++i) {
		v(i) = rand();
	}

	if (allocated) {
		delete ptr;
	}
	return v;
}

/* Pick the elements over the threshold. Return the new vector.
The count of picked elements is returned via the argument. */
template<typename Ty, typename uint_t, typename allocator_t>
DenseVector<Ty, uint_t, allocator_t> pick_threshold(
	DenseVector<Ty, uint_t, allocator_t> &v, Ty threshold, size_t& counter) {
	DenseVector<Ty, uint_t, allocator_t> v2(v.size);
	assert(threshold > 0, "Bad input.");
	counter = 0;

	for (DenseVector<Ty, uint_t, allocator_t>::size_type i = 0; i < v.size; ++i) {
		Ty elem = v(i);
		if (abs(elem) > threshold) {
			counter++;
			v2(i) = elem;
		}
	}
	return v2;
}

ns_end