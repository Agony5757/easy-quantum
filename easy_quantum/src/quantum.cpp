#include <quantum.h>
#include <omp.h>
#include <iomanip>

ns_easyquantum

//Circuit & Circuit::operator-(std::vector<fp_t> groupdef) {
//	groupdefs.push_back(groupdef);
//	return *this;
//}

void get_damping_kraus_op(
	double *k0_real, double *k0_imag,
	double *k1_real, double *k1_imag,
	const int T1, const int T_gate) {

	double ratio = T_gate * 1.0 / T1;
	double p = 1 - exp(-ratio);

	k0_real[0] = 1;
	k0_real[1] = 0;
	k0_real[2] = 0;
	k0_real[3] = sqrt(1 - p);

	k0_imag[0] = 0;
	k0_imag[1] = 0;
	k0_imag[2] = 0;
	k0_imag[3] = 0;

	k1_real[0] = 0;
	k1_real[1] = sqrt(p);
	k1_real[2] = 0;
	k1_real[3] = 0;

	k1_imag[0] = 0;
	k1_imag[1] = 0;
	k1_imag[2] = 0;
	k1_imag[3] = 0;
}

void get_dephasing_kraus_op(
	double *k0_real, double *k0_imag,
	double *k1_real, double *k1_imag,
	const int T1, const int T2, const int T_gate) {

	double p =
		(1 -
			exp(-(
				T_gate*1.0 / T2 - T_gate * 1.0 / 2 / T1
				)
			)
			) / 2;

	k0_real[0] = sqrt(1 - p);
	k0_real[1] = 0;
	k0_real[2] = 0;
	k0_real[3] = sqrt(1 - p);

	k0_imag[0] = 0;
	k0_imag[1] = 0;
	k0_imag[2] = 0;
	k0_imag[3] = 0;

	k1_real[0] = sqrt(p);
	k1_real[1] = 0;
	k1_real[2] = 0;
	k1_real[3] = -sqrt(p);

	k1_imag[0] = 0;
	k1_imag[1] = 0;
	k1_imag[2] = 0;
	k1_imag[3] = 0;
}

ns_end