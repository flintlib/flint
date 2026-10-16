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
#include "fmpq_types.h"
#include "fmpz_poly.h"
#include "longlong.h"

void revbin1(fmpz * out, const fmpz * in, slong len, slong bits);
void revbin2(fmpz * out, const fmpz * in, slong len, slong bits);
void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, slong bits);
double _fmpz_poly_evaluate_horner_d_2exp2_precomp(slong * exp, const double * poly, const slong * poly_exp, slong n, double d, slong dexp);
void _fmpz_poly_gcd_modular_primes(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, ulong first_prime);
int _checked_nmod_poly_interpolate(nn_ptr r, nn_srcptr x, nn_srcptr y, slong n, nmod_t mod);

slong _fmpz_poly_scale_positive_roots_0_1(fmpz * pol, slong len);

void _fmpz_poly_isolate_real_roots_vca(fmpq * exact_roots, slong * n_exact,
    fmpz * c_array, slong * k_array, slong * n_interval,
    const fmpz_poly_t pol, int positive_only);

int _fmpz_poly_isolate_real_roots_sturm(fmpq * exact_roots, slong * n_exact,
    fmpz * c_array, slong * k_array, slong * n_interval,
    const fmpz * pol, slong len, int positive_only, slong max_size);

int _fmpz_poly_isolate_real_roots_signs(fmpq * exact_roots, slong * n_exact,
    fmpz * c_array, slong * k_array, slong * n_interval,
    const fmpz * pol, slong len, int positive_only, slong budget);

/* Minimum length for trying root isolation using Sturm sequences
   and sign changes, respectively. */
#define FMPZ_POLY_ISOLATE_REAL_ROOTS_STURM_MIN_LEN 64
#define FMPZ_POLY_ISOLATE_REAL_ROOTS_SIGNS_MIN_LEN 32

/* Tuning parameters for real root counting ********************************/

/* Use the subresultant Sturm sequence directly when len <= MAX_LEN and
   len * bits <= MAX_SIZE. */
#define NUM_REAL_ROOTS_STURM_MAX_LEN 12
#define NUM_REAL_ROOTS_STURM_MAX_SIZE 800

/* Maximum number of evaluations for the sign change method
   (_fmpz_poly_isolate_real_roots_signs) in fmpz_poly_isolate_real_roots
   before falling back on VCA. Successful runs on real-rooted input
   typically need at most 2 len log2(len) evaluations, and about len^2 / 4
   for Eulerian polynomials (the worst case we know of); inputs with
   clusters of roots (for example Taylor shifts of Eulerian polynomials)
   make the method fail, so the budget should not be much larger. */
static inline slong
_fmpz_poly_isolate_real_roots_signs_budget(slong len)
{
    double b = FLINT_MAX(4.0 * (double) len * (double) FLINT_BIT_COUNT(len),
        (double) len * (double) len / 2.0);
    return (b >= (double) (WORD_MAX / 2)) ? WORD_MAX / 2 : (slong) b;
}

/* Give up on the primitive Sturm sequence (and use VCA) when a remainder
   has len * bits exceeding this. */
static inline slong
_fmpz_poly_num_real_roots_sturm_bound(slong len, slong bits)
{
    double b = (double) len * (double) (bits + (slong) FLINT_BIT_COUNT(len) + 8);
    return (b >= (double) (WORD_MAX / 2)) ? WORD_MAX / 2 : (slong) b;
}

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
