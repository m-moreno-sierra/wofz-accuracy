//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      hp_fref.h
//
//  Purpose:   Declares functions in namespace HP that compute the target function f_ref(z)
//             and its Taylor coefficients with arbitrary precision.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#ifndef CERF_FLINT_HP_FREF_H
#define CERF_FLINT_HP_FREF_H

#include <flint/acb.h>
#include <vector>

//! The reference function in arbitrary high precision, using data types from FLINT/Arb.

namespace HP {

//! Computes 'F' := f_ref('Z') with a precision of 'prec' decimal digits.
void acb_fref(acb_t F, const acb_t Z, slong prec);

//! Computes list of Taylor coefficient 'V' at 'Z' with a precision of 'prec' decimal digits.
void acb_fcoeffs(std::vector<acb_t>& V, const acb_t Z, slong prec);

} // namespace HP

#endif // CERF_FLINT_HP_FREF_H
