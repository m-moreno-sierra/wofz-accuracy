//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      terms.cpp
//
//  Purpose:   Provides function based on the high-precision reference function.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>
#include "terms.h"
#include "ref.h"

namespace {

double taylor_remainder(double tau, int N0, int N1, const std::vector<Ref::Coeff>& WN)
{
    assert(N1 <= WN.size());
    assert(tau > 0);
    double ret = 0;
    double tn = pow(tau, N0);
    for (int n=N0; n<N1; ++n) {
	ret += WN[n].abs * tn;
	tn *= tau;
    }
    return ret;
}

} // namespace


double Terms::truncation_error(std::complex<double> z, int N, double tau, const std::vector<Ref::Coeff>& WN)
{
    if (tau==0)
	return 0;
    static const double eps = pow(2, -53);
    const double te = ::taylor_remainder(tau, N, N+10, WN);
    const double tee = ::taylor_remainder(tau, N+10, N+20, WN);
    if (te > 1e5*eps)
	return std::numeric_limits<double>::infinity();
    if (tee > 1e-3*te) {
	std::cerr << "z = " << z << std::endl;
	std::cerr << "w[N] = " << WN[N].fn << std::endl;
	std::cerr << "te = " << te << std::endl;
	std::cerr << "tee = " << tee << std::endl;
	throw std::runtime_error("Truncation error of truncation error too large");
    }
    return te / eps;
}

double Terms::rounding_error(std::complex<double> z, int N, double tau, const std::vector<Ref::Coeff>& WN)
{
    const double lambda1 = z.real() ? sqrt(5)+1 : 2.; // lambda + 1
    double re = 0;
    double tn = 1;
    for (int n=0; n<N; ++n) {
	const int cn = (n==0 || n==N-1) ? 1 : 2;
	re += WN[n].abs * tn * (WN[n].rel_rounding_err + cn + n*lambda1);
	tn *= tau;
    }
    return re;
}
