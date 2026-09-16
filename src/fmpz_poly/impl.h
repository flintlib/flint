/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_IMPL_H
#define FMPZ_POLY_IMPL_H

#include "fmpz_types.h"

void revbin1(fmpz * out, const fmpz * in, slong len, slong bits);
void revbin2(fmpz * out, const fmpz * in, slong len, slong bits);
void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, slong bits);
double _fmpz_poly_evaluate_horner_d_2exp2_precomp(slong * exp, const double * poly, const slong * poly_exp, slong n, double d, slong dexp);
int _checked_nmod_poly_interpolate(nn_ptr r, nn_srcptr x, nn_srcptr y, slong n, nmod_t mod);

/* Internal helpers for the mpn-based multiplication routines ****************/

/* Sets {r, n len} to the coefficients of x as n-limb two's complement
   integers, adding 2^bias_bits mod 2^(FLINT_BITS n) if bias_bits >= 0. */
void _fmpz_vec_get_limbs(nn_ptr r, const fmpz * x, slong len, slong n, slong bias_bits);

/* Sets {r, (n + 1) len} to the coefficients of x in sign-magnitude form
   (a sign limb followed by the n-limb magnitude). */
void _fmpz_vec_get_limbs_signmag(nn_ptr r, const fmpz * x, slong len, slong n);

/* Given the product coefficients U[i] for nlo <= i < nhi (of the packed
   polynomials a and b with n1 and n2 limbs), stored contiguously with
   stride s starting from index nlo, sets res[i - nlo] to the correct
   fmpz values, applying the bias correction if needed. The method is one
   of the FLINT_MPN_POLY_MUL_* representations. */
void _fmpz_poly_mpn_set_outputs(fmpz * res, nn_ptr U, slong s, slong nlo, slong nhi, int method, nn_srcptr a, slong len1, slong n1, nn_srcptr b, slong len2, slong n2, slong ub1, slong ub2, int sgn1, int sgn2, int squaring);

#endif
