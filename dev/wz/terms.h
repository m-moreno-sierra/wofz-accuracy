//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      terms.h
//
//  Purpose:   Declares function based on the high-precision reference function.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#ifndef CERF_FLINT_TERMS_H
#define CERF_FLINT_TERMS_H

#include "wn.h"

namespace Terms {

double wmin(int ix, int iy, int d2, const std::vector<std::vector<double>>& W);

double truncation_error(std::complex<double> z, int N, double tau,
			const std::vector<Ref::Coeff>& WN);

double rounding_error(std::complex<double> z, int N, double tau,
		      const std::vector<Ref::Coeff>& WN);

} // namespace Terms

#endif // CERF_FLINT_TERMS_H
