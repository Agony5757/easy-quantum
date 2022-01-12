#pragma once
#include "mathutil.h"
#include "global.h"
#include "randomutil.h"

ns_easyquantum

template<typename qtraits_t = default_qtraits>
struct state_manipulator {
	using qtraits = qtraits_t;
	using fp_t = typename qtraits::value_t;
	using qid = typename qtraits::qidx_t;
	using uint_t = typename qtraits::idx_t;

	static uint_t create_state(fp_t** real, fp_t** imag, const qid qn) {
		*real = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		*imag = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		memset(*real, 0, (sizeof(fp_t)) * (pow2(qn)));
		memset(*imag, 0, (sizeof(fp_t)) * (pow2(qn)));
		*real[0] = 1;
		return pow2(qn);
	}

	static void copy_state(fp_t** real_to, fp_t** imag_to,
		fp_t* real_from, fp_t* imag_from, const uint_t size) {

		*real_to = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		*imag_to = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));

		for (uint_t i = 0; i < size; ++i) {
			*real_to[i] = real_from[i];
			*imag_to[i] = imag_from[i];
		}
	}

	static void refresh_state(fp_t* real, fp_t* imag, const uint_t size) {
		memset(real, 0, sizeof(fp_t) * size);
		memset(imag, 0, sizeof(fp_t) * size);
		real[0] = 1;
	}

	static uint_t create_state_full_hadamard(fp_t** real, fp_t** imag, const qid qn)
	{
		*real = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		*imag = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		fp_t p = sqrt(1.0 / pow2(qn));
		for (uint_t i = 0; i < pow2(qn); ++i) {
			(*real)[i] = p;
			(*imag)[i] = 0;
		}
		return pow2(qn);
	}

	static uint_t create_state_randomized(fp_t** real, fp_t** imag, const qid qn)
	{
		*real = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));
		*imag = (fp_t*)malloc(sizeof(fp_t) * (pow2(qn)));

		for (uint_t i = 0; i < pow2(qn); ++i) {
			(*real)[i] = (fp_t)_default_random_generator();
			(*imag)[i] = (fp_t)_default_random_generator();
		}

		// normalize
		fp_t abs = 0.0;
		for (uint_t i = 0; i < pow2(qn); ++i) {
			abs += (*real)[i] * (*real)[i];
			abs += (*imag)[i] * (*imag)[i];
		}

		fp_t norm2 = sqrt(abs);
		for (uint_t i = 0; i < pow2(qn); ++i) {
			(*real)[i] /= norm2;
			(*imag)[i] /= norm2;
		}

		return pow2(qn);
	}

	static void output_state(const fp_t* real, const fp_t* imag, const uint_t size) {
		for (uint_t i = 0; i < size; ++i) {
			printf("state[%lld] = %lf + %lf i\n", i, real[i], imag[i]);
		}
		printf("\n");
	}

	static void state_to_file(const fp_t* real, const fp_t* imag, const uint_t size,
		std::ostream& out) {
		for (uint_t i = 0; i < size; ++i) {
			char s[100] = { '\0' };
			sprintf(s, "state[%lld] = %lf + %lf i\n", i, real[i], imag[i]);
			out << s;
		}
		out << std::endl;
	}

	static void print_gate1q(const fp_t* real, const fp_t* imag) {
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				int idx = i * 2 + j;
				if (round0(imag[idx], quantum_default_threshold))
					printf("%lf ", real[idx]);
				else if (round0(real[idx], quantum_default_threshold))
					printf("%lfi ", imag[idx]);
				else
					printf("%lf+%lfi ", real[idx], imag[idx]);
			}
			printf("\n");
		}
		printf("\n");
	}

	static fp_t GetNorm2i(const fp_t* real, const fp_t* imag, uint_t i) {
		return real[i] * real[i] + imag[i] * imag[i];
	}

	static int u1q_mul_vec(fp_t* real, fp_t* imag, const uint_t size,
		const fp_t* u_real, const fp_t* u_imag, const qid qn) {

		fp_t u00r = u_real[0], u01r = u_real[1], u10r = u_real[2], u11r = u_real[3];
		fp_t u00i = u_imag[0], u01i = u_imag[1], u10i = u_imag[2], u11i = u_imag[3];

		for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
			for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
				fp_t r0 = real[j], i0 = imag[j],
					r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
				fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
					res_temp0i, res_temp1i, res_temp2i, res_temp3i;

				complex_multiplication(res_temp0r, res_temp0i, u00r, u00i, r0, i0);
				complex_multiplication(res_temp1r, res_temp1i, u01r, u01i, r1, i1);
				complex_multiplication(res_temp2r, res_temp2i, u10r, u10i, r0, i0);
				complex_multiplication(res_temp3r, res_temp3i, u11r, u11i, r1, i1);
				r0 = res_temp0r + res_temp1r; real[j] = r0;
				i0 = res_temp0i + res_temp1i; imag[j] = i0;
				r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1;
				i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1;
			}
		}
		return 0;
	}

	static int cu_mul_vec(fp_t* real, fp_t* imag, const uint_t size,
		const fp_t* u_real, const fp_t* u_imag, const qid control, const qid qn) {
		fp_t u00r = u_real[0], u01r = u_real[1], u10r = u_real[2], u11r = u_real[3];
		fp_t u00i = u_imag[0], u01i = u_imag[1], u10i = u_imag[2], u11i = u_imag[3];

		for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
			for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
				if (!((j >> control) % 2)) continue;
				fp_t r0 = real[j], i0 = imag[j],
					r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
				fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
					res_temp0i, res_temp1i, res_temp2i, res_temp3i;

				complex_multiplication(res_temp0r, res_temp0i, u00r, u00i, r0, i0);
				complex_multiplication(res_temp1r, res_temp1i, u01r, u01i, r1, i1);
				complex_multiplication(res_temp2r, res_temp2i, u10r, u10i, r0, i0);
				complex_multiplication(res_temp3r, res_temp3i, u11r, u11i, r1, i1);
				r0 = res_temp0r + res_temp1r; real[j] = r0;
				i0 = res_temp0i + res_temp1i; imag[j] = i0;
				r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1;
				i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1;
			}
		}
		return 0;
	}

	static int cnot_mul_vec(fp_t* real, fp_t* imag, const uint_t size,
		const qid control, const qid qn) {

		for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
			for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
				if (!((j >> control) % 2)) continue;
				fp_t r0 = real[j], i0 = imag[j],
					r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
				fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
					res_temp0i, res_temp1i, res_temp2i, res_temp3i;

				complex_multiplication(res_temp0r, res_temp0i, 0.0, 0.0, r0, i0);
				complex_multiplication(res_temp1r, res_temp1i, 1.0, 0.0, r1, i1);
				complex_multiplication(res_temp2r, res_temp2i, 1.0, 0.0, r0, i0);
				complex_multiplication(res_temp3r, res_temp3i, 0.0, 0.0, r1, i1);
				r0 = res_temp0r + res_temp1r; real[j] = r0;
				i0 = res_temp0i + res_temp1i; imag[j] = i0;
				r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1;
				i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1;
			}
		}
		return 0;
	}

	static int ccnot_mul_vec(fp_t* real, fp_t* imag, const uint_t size,
		const qid control1, const qid control2, const qid qn) {

		for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
			for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
				if (!((j >> control1) % 2)) continue;
				if (!((j >> control2) % 2)) continue;
				fp_t r0 = real[j], i0 = imag[j],
					r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
				fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
					res_temp0i, res_temp1i, res_temp2i, res_temp3i;

				complex_multiplication(res_temp0r, res_temp0i, 0, 0, r0, i0);
				complex_multiplication(res_temp1r, res_temp1i, 1, 0, r1, i1);
				complex_multiplication(res_temp2r, res_temp2i, 1, 0, r0, i0);
				complex_multiplication(res_temp3r, res_temp3i, 0, 0, r1, i1);
				r0 = res_temp0r + res_temp1r; real[j] = r0;
				i0 = res_temp0i + res_temp1i; imag[j] = i0;
				r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1;
				i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1;
			}
		}
		return 0;
	}

	static int diagonal_mul_vec(fp_t* real, fp_t* imag, const uint_t size,
		const fp_t* h, const fp_t theta) {
		for (uint_t i = 0; i < size; ++i) {
			fp_t res_real;
			fp_t res_imag;
			fp_t real_angle = cos(h[i] * theta);
			fp_t imag_angle = -sin(h[i] * theta);
			complex_multiplication(res_real, res_imag, real_angle, imag_angle, real[i], imag[i]);
			real[i] = res_real;
			imag[i] = res_imag;
		}
		return 0;
	}

	static int vec_boundary_prob(fp_t* prob, const fp_t* real, const fp_t* imag,
		const uint_t size, const qid* qn, const qid qn_size) {
		uint_t prob_size = pow2(qn_size);

		for (uint_t i = 0; i < size; ++i) {
			// 对于任何一个i，使之按照qn进行排列
			uint_t x = 0;
			for (qid j = 0; j < qn_size; ++j) {
				x += (((i >> qn[j]) % 2) << j);
			}
			prob[x] += GetNorm2i(real, imag, i);
		}
		return 0;
	}

	static uint_t pick_digit(const uint_t index, const qid* qn, const qid qn_size) {
		uint_t new_index = 0;
		for (qid i = 0; i < qn_size; ++i) {
			qid bit = qn[i];
			bool one = (index >> bit) & 1;
			new_index += (one ? pow2(i) : 1);
		}
		return new_index;
	}

	static fp_t* get_prob(const fp_t* real, const fp_t* imag, const uint_t size) {

		fp_t* probs = (fp_t*)malloc(sizeof(fp_t) * size);
		// memset(probs, 0, sizeof(fp_t) * size);
		for (uint_t i = 0; i < size; ++i) {
			probs[i] = GetNorm2i(real, imag, i);
		}
		return probs;
	}

	static void get_prob(fp_t* probs, const fp_t* real, const fp_t* imag, const uint_t size) {
		for (uint_t i = 0; i < size; ++i) {
			probs[i] = GetNorm2i(real, imag, i);
		}
	}

	static fp_t* get_partial_prob(const fp_t* real, const fp_t* imag, const uint_t size, const qid* qn, const qid qn_size) {
		uint_t probs_size = pow2(qn_size);
		fp_t* probs = (fp_t*)malloc(sizeof(fp_t) * probs_size);
		memset(probs, 0, sizeof(fp_t) * probs_size);
		for (uint_t i = 0; i < size; ++i) {
			uint_t picked = pick_digit(i, qn, qn_size);
			probs[picked] += (GetNorm2i(real, imag, i));
		}
		return probs;
	}

	static void get_partial_prob(fp_t* probs, const fp_t* real, const fp_t* imag, const uint_t size, const qid* qn, const qid qn_size) {
		for (uint_t i = 0; i < size; ++i) {
			uint_t picked = pick_digit(i, qn, qn_size);
			probs[picked] += (GetNorm2i(real, imag, i));
		}
	}

	static uint_t measure(fp_t* real, fp_t* imag, const uint_t size, const qid* qn, const qid qn_size, RandomEngine* rng) {
		if (qn_size == 1) {
			return measure_1q(real, imag, size, qn[0], rng);
		}
		else if (qn_size > 1) {
			fp_t* probs = get_partial_prob(real, imag, size, qn, qn_size);
			fp_t random = (*rng)();
			fp_t accumulate = 0.;

			uint_t meas_res = 0;

			for (uint_t i = 0; i < size; ++i) {
				accumulate += probs[i];
				if (random < accumulate) {
					meas_res = i;
					break;
				}
			}
			for (uint_t i = 0; i < size; ++i) {
				if (pick_digit(i, qn, qn_size) != meas_res) {
					real[i] = 0;
					imag[i] = 0;
				}
				else {
					real[i] /= probs[meas_res];
					imag[i] /= probs[meas_res];
				}
			}
			free(probs);
			return meas_res;
		}
		else {
			return -1;
		}
	}

	static uint_t const_measure(const fp_t* real, const fp_t* imag,
		const uint_t size, const qid* qn, const qid qn_size, RandomEngine* rng) {
		// Non-demolition measurement
		if (qn_size == 1) {
			return const_measure_1q(real, imag, size, qn[0], rng);
		}
		else if (qn_size > 1) {
			fp_t* probs = get_partial_prob(real, imag, size, qn, qn_size);
			fp_t random = (*rng)();
			fp_t accumulate = 0.;
			uint_t meas_res = 0;

			for (uint_t i = 0; i < size; ++i) {
				accumulate += probs[i];
				if (random < accumulate) {
					meas_res = i;
					break;
				}
			}
			free(probs);
			return meas_res;
		}
		else {
			return -1;
		}
	}

	static uint_t measure_1q(fp_t* real, fp_t* imag, const uint_t size,
		const qid qn, RandomEngine* rng) {
		fp_t prob0 = 0;
		uint_t res = 0;

		for (uint_t i = 0; i < size; ++i) {
			if ((size >> qn) % 2 == 0) {
				prob0 += (GetNorm2i(real, imag, i));
			}
		}
		fp_t random = (*rng)();
		if (random < prob0) {
			for (uint_t i = 0; i < size; ++i) {
				if ((size >> qn) % 2 == 0) {
					real[i] /= prob0;
					imag[i] /= prob0;
				}
			}
			res = 0;
		}
		else {
			fp_t prob1 = 1 - prob0;
			for (uint_t i = 0; i < size; ++i) {
				if ((size >> qn) % 2 == 0) {
					real[i] /= prob1;
					imag[i] /= prob1;
				}
			}
			res = 1;
		}
		return res;
	}

	static uint_t const_measure_1q(const fp_t* real, const fp_t* imag,
		const uint_t size, const qid qn, RandomEngine* rng) {
		fp_t prob0 = 0;
		uint_t res = 0;

		for (uint_t i = 0; i < size; ++i) {
			if ((size >> qn) % 2 == 0) {
				prob0 += (GetNorm2i(real, imag, i));
			}
		}
		fp_t random = (*rng)();
		if (random < prob0) {
			res = 0;
		}
		else {
			res = 1;
		}
		return res;
	}

	static uint_t binary_search(const fp_t* arr, const uint_t size, const fp_t r) {
		// suppose arr is sorted from small to big

		uint_t range_left = 0;
		uint_t interval = size / 2;
		uint_t range_right = size;
		while (1) {
			if (arr[interval] == r)
				return interval;
			else if (arr[interval] > r) {
				// in left half
				uint_t new_interval = (interval + range_left) / 2;
				if (new_interval == range_left) {
					return new_interval;
				}
				range_right = interval;
				interval = new_interval;
			}
			else if (arr[interval] < r) {
				uint_t new_interval = (interval + range_right) / 2;
				if (new_interval == interval) {
					return new_interval;
				}
				range_left = interval;
				interval = new_interval;
			}
		}
	}

	static uint_t* measure_n_shots_thread(
		const fp_t* real, const fp_t* imag, const uint_t size,
		const qid* qn, const qid qn_size, RandomEngine* rng, uint_t shots) {

		if (size == pow2(qn_size)) {
			return measure_n_shots_all_thread(real, imag, size, rng, shots);
		}
		else if (qn_size > 1) {
			// first calculate the size of accumulated probs
			uint_t accum_prob_size = (size << qn_size);
			assert(accum_prob_size > 0, "Bad size.");

			fp_t* accum_prob = (fp_t*)malloc(sizeof(fp_t) * accum_prob_size);
			uint_t res_size = pow2(qn_size);
			uint_t* res = (uint_t*)malloc(sizeof(uint_t) * shots);
			accum_prob[0] = 0;
			fp_t* probs = get_partial_prob(real, imag, size, qn, qn_size);
			for (uint_t i = 0; i < accum_prob_size - 1; ++i) {
				accum_prob[i + 1] = accum_prob[i] + probs[i];
			}
#pragma omp parallel
			for (int i = 0; i < shots; ++i) {
				fp_t random = (*rng)();
				res[i] = binary_search(accum_prob, size, random);
			}
			free(probs);
			free(accum_prob);
			return res;
		}
		else {
			assert(false, "Bad size.");
		}
	}

	static uint_t* measure_n_shots_all_thread(const fp_t* real, const fp_t* imag,
		const uint_t size, RandomEngine* rng, uint_t shots) {

		// first get accumulated probs
		// then use binary search
		fp_t* accum_prob = (fp_t*)malloc(sizeof(fp_t) * size);
		uint_t* res = (uint_t*)malloc(sizeof(uint_t) * shots);
		accum_prob[0] = 0;
		for (uint_t i = 0; i < size - 1; ++i) {
			accum_prob[i + 1] = accum_prob[i] + GetNorm2i(real, imag, i);
		}
#pragma omp parallel
		for (int i = 0; i < shots; ++i) {
			fp_t random = (*rng)();
			res[i] = binary_search(accum_prob, size, random);
		}

		free(accum_prob);
		return res;
	}

	static int rx(fp_t* real, fp_t* imag, const uint_t size, const fp_t theta, const qid qn) {
		fp_t u_real[4] = { cos(theta / 2), 0, 0, cos(theta / 2) };
		fp_t u_imag[4] = { 0, -sin(theta / 2), -sin(theta / 2), 0 };

		return u1q_mul_vec(real, imag, size, u_real, u_imag, qn);
	}

	static int ry(fp_t* real, fp_t* imag, const uint_t size, const fp_t theta, const qid qn) {
		fp_t u_real[4] = { cos(theta / 2), -sin(theta / 2), sin(theta / 2), cos(theta / 2) };
		fp_t u_imag[4] = { 0, 0, 0, 0 };

		return u1q_mul_vec(real, imag, size, u_real, u_imag, qn);
	}

	static int rz(fp_t* real, fp_t* imag, const uint_t size, const fp_t theta, const qid qn) {
		fp_t u_real[4] = { 1, 0, 0, cos(theta) };
		fp_t u_imag[4] = { 0, 0, 0, sin(theta) };

		return u1q_mul_vec(real, imag, size, u_real, u_imag, qn);
	}

	static int h(fp_t* real, fp_t* imag, const uint_t size, const qid qn) {
		fp_t u_real[4] = { 1 / SQ2, 1 / SQ2, 1 / SQ2, -1 / SQ2 };
		fp_t u_imag[4] = { 0, 0, 0, 0 };

		return u1q_mul_vec(real, imag, size, u_real, u_imag, qn);
	}

	static int cnot(fp_t* real, fp_t* imag, const uint_t size, const qid control, const qid target) {
		return cnot_mul_vec(real, imag, size, control, target);
	}

	static int perform_kraus(fp_t* real, fp_t* imag, const uint_t size,
		const double* k0real, const double* k0imag,
		const double* k1real, const double* k1imag,
		const qid qn, const fp_t random) {

		/* get p0 */
		double p0 = 0;

		fp_t u00r = k0real[0], u01r = k0real[1], u10r = k0real[2], u11r = k0real[3];
		fp_t u00i = k0imag[0], u01i = k0imag[1], u10i = k0imag[2], u11i = k0imag[3];

		for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
			for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
				fp_t r0 = real[j], i0 = imag[j],
					r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
				fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
					res_temp0i, res_temp1i, res_temp2i, res_temp3i;

				complex_multiplication(res_temp0r, res_temp0i, u00r, u00i, r0, i0);
				complex_multiplication(res_temp1r, res_temp1i, u01r, u01i, r1, i1);
				complex_multiplication(res_temp2r, res_temp2i, u10r, u10i, r0, i0);
				complex_multiplication(res_temp3r, res_temp3i, u11r, u11i, r1, i1);
				r0 = res_temp0r + res_temp1r;
				i0 = res_temp0i + res_temp1i;
				r1 = res_temp2r + res_temp3r;
				i1 = res_temp2i + res_temp3i;

				p0 += (r0 * r0 + i0 * i0 + r1 * r1 + i1 * i1);
			}
		}
		//debug_info_s(Kraus Performing : );
		//debug_info_s(KrausOp1);
		//debug_print_mat22(k0real, k0imag);
		//debug_info_s(KrausOp2);
		//debug_print_mat22(k1real, k1imag);
		//debug_display(p0);
		//debug_display(random);

		if (random < p0) {
			for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
				for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
					fp_t r0 = real[j], i0 = imag[j],
						r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
					fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
						res_temp0i, res_temp1i, res_temp2i, res_temp3i;

					complex_multiplication(res_temp0r, res_temp0i, u00r, u00i, r0, i0);
					complex_multiplication(res_temp1r, res_temp1i, u01r, u01i, r1, i1);
					complex_multiplication(res_temp2r, res_temp2i, u10r, u10i, r0, i0);
					complex_multiplication(res_temp3r, res_temp3i, u11r, u11i, r1, i1);
					r0 = res_temp0r + res_temp1r; real[j] = r0 / sqrt(p0);
					i0 = res_temp0i + res_temp1i; imag[j] = i0 / sqrt(p0);
					r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1 / sqrt(p0);
					i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1 / sqrt(p0);
				}
			}
		}
		else {
			fp_t u00r = k1real[0], u01r = k1real[1], u10r = k1real[2], u11r = k1real[3];
			fp_t u00i = k1imag[0], u01i = k1imag[1], u10i = k1imag[2], u11i = k1imag[3];

			for (uint_t i = 0; i < size; i += (pow2(qn + 1))) {
				for (uint_t j = i; j < (i + (pow2(qn))); ++j) {
					fp_t r0 = real[j], i0 = imag[j],
						r1 = real[j + (pow2(qn))], i1 = imag[j + (pow2(qn))];
					fp_t res_temp0r, res_temp1r, res_temp2r, res_temp3r,
						res_temp0i, res_temp1i, res_temp2i, res_temp3i;

					complex_multiplication(res_temp0r, res_temp0i, u00r, u00i, r0, i0);
					complex_multiplication(res_temp1r, res_temp1i, u01r, u01i, r1, i1);
					complex_multiplication(res_temp2r, res_temp2i, u10r, u10i, r0, i0);
					complex_multiplication(res_temp3r, res_temp3i, u11r, u11i, r1, i1);
					r0 = res_temp0r + res_temp1r; real[j] = r0 / sqrt(1 - p0);
					i0 = res_temp0i + res_temp1i; imag[j] = i0 / sqrt(1 - p0);
					r1 = res_temp2r + res_temp3r; real[j + (pow2(qn))] = r1 / sqrt(1 - p0);
					i1 = res_temp2i + res_temp3i; imag[j + (pow2(qn))] = i1 / sqrt(1 - p0);
				}
			}
		}
		//debug_output_state;
		return 0;
	}
};

template<typename qtraits_t = default_qtraits>
// typename value_t, typename idx_t, typename qidx_t, qidx_t qn>
struct state {
	using qtraits = qtraits_t;
	using qid = typename qtraits::qidx_t;
	using value_t = typename qtraits::value_t;
	using uint_t = typename qtraits::idx_t;

	value_t* real = nullptr;
	value_t* imag = nullptr;
	uint_t size = 0;

	state() = delete;
	state(qid qn) {
		size = state_manipulator<qtraits>::create_state(&real, &imag, qn);
	}

	void refresh() {
		state_manipulator<qtraits>::refresh_state(real, imag, size);
	}
};

ns_end