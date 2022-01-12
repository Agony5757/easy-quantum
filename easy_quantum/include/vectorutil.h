#pragma once
#include <global.h>
#include <testutil.h>

ns_easyquantum

/* Check if an elem is in a vector, return boolean. */
template<typename T>
bool is_in(const std::vector<T> &vec, const T &elem) {
	if (find(vec.begin(), vec.end(), elem) != vec.end()) {
		return true;
	}
	else return false;
}

/* Check if an elem is in a vector, return boolean. */
template<typename T, typename T2, typename Func>
bool is_in(const std::vector<T>& vec, const T2& elem, Func& f) {
	for (const T& e : vec){
		// whenever find return true
		if (f(e, elem)) return true;
	}
	// cannot find
	return false;
}

/* Erase elements in another vector, return a new vector */
template<typename T>
std::vector<T> erase_many(std::vector<T> vec, const std::vector<T> &to_erase) {
	std::vector<T> vec1(vec.begin(), vec.end());
	auto im = std::remove_if(vec1.begin(), vec1.end(),
		[&to_erase](T& elem) {return is_in(to_erase, elem); });
	
	vec1.erase(im, vec1.end());
	return vec1;
}

/* For all elements do function, based on for_each */
template<class _InIt, class _Fn> 
_Fn for_all(_InIt &_Capacitor, _Fn _Func) {
	return for_each(_Capacitor.begin(), _Capacitor.end(), _Func);
}

/* Check if two containers are identical, return boolean. 
Suppose inputs are containers. */
template<typename T>
bool is_same_container(const T& a, const T& b) {
	if (a.size() != b.size()) {
		return false;
	}
	for (size_t i = 0; i < a.size(); ++i) {
		if (a[i] != b[i])
			return false;
	}
	return true;
}

/* Append one vector with another. */
template<typename T>
void merge_into(std::vector<T> &to, const std::vector<T> &from) {
	to.insert(to.end(), from.begin(), from.end());
}

/* Type tools. Obtain a type's vector wrapper. */
template<typename Ty>
constexpr std::vector<Ty> _to_vec(Ty) { return std::vector<Ty>(); }

template<typename Ty>
struct is_numerical {
	static constexpr bool value = std::is_integral_v<Ty> | std::is_floating_point_v<Ty>;
};

template<typename Ty> constexpr bool is_numerical_v = is_numerical<Ty>::value;

/* X' = aX */
template<typename Ty, typename ConstTy>
auto VecMulConst(std::vector<Ty> vec, ConstTy constant) -> decltype(_to_vec(Ty()*ConstTy())) {
	static_assert(is_numerical_v<Ty>, "Bad Type.");
	static_assert(is_numerical_v<ConstTy>, "Bad Type.");

	std::vector<decltype(Ty()*ConstTy())> res_vec(vec.size());
	ForAlli(vec) { res_vec[i] = vec[i] * constant; }
	return res_vec;
}

/* X' = aX */
template<typename Ty, typename ConstTy>
auto operator*(std::vector<Ty> vec, ConstTy constant) -> decltype(VecMulConst(vec, constant)) {
	return VecMulConst<Ty, ConstTy>(vec, constant);
}

/* X' = aX */
template<typename Ty, typename ConstTy>
auto operator*(ConstTy constant, std::vector<Ty> vec) -> decltype(VecMulConst(vec, constant)) {
	return VecMulConst<Ty, ConstTy>(vec, constant);
}

/* X3 = X1.*X2 */
template<typename Ty1, typename Ty2>
auto VecInnerProduct(std::vector<Ty1> vec1, std::vector<Ty2> vec2)
-> decltype(_to_vec(Ty1()*Ty2())) {

	static_assert(is_numerical_v<Ty1>, "Bad Type.");
	static_assert(is_numerical_v<Ty2>, "Bad Type.");

	std::vector<decltype(Ty1()*Ty2())> res_vec(vec1.size());
	assert(vec1.size() == vec2.size(), "Bad Type.");
	if (vec1.size() == vec2.size()) {
		ForAlli(vec1) {
			res_vec[i] += (vec1[i] * vec2[i]);
		}
		return res_vec;
	}
}

/* d = SUM(X) */
template<typename Ty>
Ty VecSum(std::vector<Ty> vec) {
	static_assert(is_numerical_v<Ty>, "Bad Type.");

	Ty data = 0;
	ForAlli(vec) {
		data += vec[i];
	}
	return data;
}

/* d = norm2(X) */
template<typename Ty>
Ty VecNorm2(std::vector<Ty> vec) {
	return VecSum(VecInnerProduct<Ty, Ty>(vec, vec));
}

/* X' = a+X */
template<typename Ty, typename ConstTy>
auto VecAddConst(std::vector<Ty> vec, ConstTy constant) -> decltype(_to_vec(Ty() + ConstTy())) {
	static_assert(is_numerical_v<Ty>, "Bad Type.");
	static_assert(is_numerical_v<ConstTy>, "Bad Type.");
	std::vector<decltype(Ty1() + Ty2())> res_vec(vec.size());
	ForAlli(vec) { res_vec[i] = vec[i] + constant; }
	return res_vec;
}

/* X' = a+X */
template<typename Ty, typename ConstTy>
auto operator+(std::vector<Ty> vec, ConstTy constant) -> decltype(VecAddConst(vec, constant)) {
	return VecAddConst<Ty, ConstTy>(vec, constant);
}

/* X' = a+X */
template<typename Ty, typename ConstTy>
auto operator+(ConstTy constant, std::vector<Ty> vec) -> decltype(VecAddConst(vec, constant)) {
	return VecAddConst<Ty, ConstTy>(vec, constant);
}

/* X' = X+Y */
template<typename Ty1, typename Ty2>
auto VecAddVec(std::vector<Ty1> vec1, std::vector<Ty2> vec2) -> decltype(_to_vec(Ty1() + Ty2())) {
	static_assert(is_numerical_v<Ty1>, "Bad Type.");
	static_assert(is_numerical_v<Ty2>, "Bad Type.");
	assert(vec1.size() == vec2.size(), "Bad Size.");

	std::vector<decltype(Ty1() + Ty2())> res_vec(vec1.size());
	ForAlli(res_vec) { res_vec[i] = vec1[i] + vec2[i]; }
	return res_vec;
}

template<typename Ty1, typename Ty2>
auto operator+(std::vector<Ty1> vec1, std::vector<Ty2> vec2) -> decltype(_to_vec(Ty1() + Ty2())) {
	return VecAddVec<Ty1, Ty2>(vec1, vec2);
}

/* Accumulate over a vector */
template<typename Ty>
Ty accumulative_addition(const std::vector<Ty> &vec) {
	Ty s = 0;
	ForAll(vec, i) {
		s += vec[i];
	}
	return s;
}

/* Accumulate over a vector with a function */
template<typename Ty>
Ty accumulative_addition(const std::vector<Ty>& vec, std::function<Ty(Ty)>& f) {
	Ty s = 0;
	ForAll(vec, i) {
		s += f(vec[i]);
	}
	return s;
}

/* Accumulative multiplication */
template<typename Ty>
Ty accumulative_multiplication(const std::vector<Ty>& vec) {
	Ty s = 1;
	ForAll(vec, i) {
		s *= vec[i];
	}
	return s;
}


ns_end