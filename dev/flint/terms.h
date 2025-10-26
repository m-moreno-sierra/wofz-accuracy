#ifndef CERF_FLINT_TERMS_H
#define CERF_FLINT_TERMS_H

#include "wn.h"

std::vector<std::vector<double>> wOnGrid(double inv_b, int d2max);
double wmin(int ix, int iy, int d2, const std::vector<std::vector<double>>& W);

double truncation_error(std::complex<double> z, int N, double tau,
			const std::vector<Coeff>& WN);

double rounding_error(std::complex<double> z, int N, double tau,
		      const std::vector<Coeff>& WN);

#endif // CERF_FLINT_TERMS_H
