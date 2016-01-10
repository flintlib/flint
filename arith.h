/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/******************************************************************************

 Copyright (C) 2010-2012 Fredrik Johansson

******************************************************************************/

#ifndef ARITH_H
#define ARITH_H

#ifdef ARITH_INLINES_C
#define ARITH_INLINE FLINT_DLL
#else
#define ARITH_INLINE static __inline__
#endif

#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"
#include "double_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* MPFR extras ***************************************************************/

FLINT_DLL void mpfr_zeta_inv_euler_product(mpfr_t res, ulong s, int char_4);

/* Various arithmetic functions **********************************************/

#define arith_primorial fmpz_primorial

FLINT_DLL void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n);
FLINT_DLL void arith_harmonic_number(fmpq_t x, slong n);

FLINT_DLL void arith_ramanujan_tau(fmpz_t res, const fmpz_t n);
FLINT_DLL void arith_ramanujan_tau_series(fmpz_poly_t res, slong n);

FLINT_DLL void arith_divisors(fmpz_poly_t res, const fmpz_t n);

#define arith_divisor_sigma fmpz_divisor_sigma
#define arith_moebius_mu fmpz_moebius_mu
#define arith_euler_phi fmpz_euler_phi

/* Stirling numbers **********************************************************/

FLINT_DLL void arith_stirling_number_1u(fmpz_t s, slong n, slong k);
FLINT_DLL void arith_stirling_number_1(fmpz_t s, slong n, slong k);
FLINT_DLL void arith_stirling_number_2(fmpz_t s, slong n, slong k);

FLINT_DLL void arith_stirling_number_1u_vec(fmpz * row, slong n, slong klen);
FLINT_DLL void arith_stirling_number_1_vec(fmpz * row, slong n, slong klen);
FLINT_DLL void arith_stirling_number_2_vec(fmpz * row, slong n, slong klen);

FLINT_DLL void arith_stirling_number_1u_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);
FLINT_DLL void arith_stirling_number_1_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);
FLINT_DLL void arith_stirling_number_2_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);

FLINT_DLL void arith_stirling_matrix_1u(fmpz_mat_t mat);
FLINT_DLL void arith_stirling_matrix_1(fmpz_mat_t mat);
FLINT_DLL void arith_stirling_matrix_2(fmpz_mat_t mat);

/* Bell numbers **************************************************************/

#if FLINT64
#define BELL_NUMBER_TAB_SIZE 26
#else
#define BELL_NUMBER_TAB_SIZE 16
#endif

extern const mp_limb_t bell_number_tab[];

FLINT_DLL double arith_bell_number_size(ulong n);

FLINT_DLL void arith_bell_number(fmpz_t b, ulong n);
FLINT_DLL void arith_bell_number_bsplit(fmpz_t res, ulong n);
FLINT_DLL void arith_bell_number_multi_mod(fmpz_t res, ulong n);

FLINT_DLL void arith_bell_number_vec(fmpz * b, slong n);
FLINT_DLL void arith_bell_number_vec_recursive(fmpz * b, slong n);
FLINT_DLL void arith_bell_number_vec_multi_mod(fmpz * b, slong n);

FLINT_DLL mp_limb_t arith_bell_number_nmod(ulong n, nmod_t mod);

FLINT_DLL void arith_bell_number_nmod_vec(mp_ptr b, slong n, nmod_t mod);
FLINT_DLL void arith_bell_number_nmod_vec_recursive(mp_ptr b, slong n, nmod_t mod);
FLINT_DLL void arith_bell_number_nmod_vec_series(mp_ptr b, slong n, nmod_t mod);


/* Euler numbers *************************************************************/

#if FLINT64
#define SMALL_EULER_LIMIT 25
#else
#define SMALL_EULER_LIMIT 15
#endif

static const mp_limb_t euler_number_small[] = {
    UWORD(1), UWORD(1), UWORD(5), UWORD(61), UWORD(1385), UWORD(50521), UWORD(2702765),
    UWORD(199360981),
#if FLINT64
    UWORD(19391512145), UWORD(2404879675441), UWORD(370371188237525),
    UWORD(69348874393137901), UWORD(15514534163557086905)
#endif
};

FLINT_DLL double arith_euler_number_size(ulong n);

FLINT_DLL void arith_euler_number_vec(fmpz * res, slong n);

FLINT_DLL void _arith_euler_number_zeta(fmpz_t res, ulong n);
FLINT_DLL void arith_euler_number(fmpz_t res, ulong n);

FLINT_DLL void arith_euler_polynomial(fmpq_poly_t poly, ulong n);

/* Bernoulli numbers *********************************************************/

#if FLINT64
#define BERNOULLI_SMALL_NUMER_LIMIT 35
#else
#define BERNOULLI_SMALL_NUMER_LIMIT 27
#endif

static const slong _bernoulli_numer_small[] = {
    WORD(1), WORD(1), WORD(-1), WORD(1), WORD(-1), WORD(5), WORD(-691), WORD(7), WORD(-3617), WORD(43867), WORD(-174611), WORD(854513),
    WORD(-236364091), WORD(8553103),
#if FLINT64
    WORD(-23749461029), WORD(8615841276005), WORD(-7709321041217), WORD(2577687858367)
#endif
};

FLINT_DLL void _arith_bernoulli_number(fmpz_t num, fmpz_t den, ulong n);
FLINT_DLL void arith_bernoulli_number(fmpq_t x, ulong n);

FLINT_DLL void _arith_bernoulli_number_vec(fmpz * num, fmpz * den, slong n);
FLINT_DLL void arith_bernoulli_number_vec(fmpq * num, slong n);

FLINT_DLL void arith_bernoulli_number_denom(fmpz_t den, ulong n);
FLINT_DLL double arith_bernoulli_number_size(ulong n);

FLINT_DLL void arith_bernoulli_polynomial(fmpq_poly_t poly, ulong n);

FLINT_DLL void _arith_bernoulli_number_zeta(fmpz_t num, fmpz_t den, ulong n);
FLINT_DLL void _arith_bernoulli_number_vec_multi_mod(fmpz * num, fmpz * den, slong n);
FLINT_DLL void _arith_bernoulli_number_vec_recursive(fmpz * num, fmpz * den, slong n);
FLINT_DLL void _arith_bernoulli_number_vec_zeta(fmpz * num, fmpz * den, slong n);

/* Cyclotomic polynomials ****************************************************/

#define _arith_cyclotomic_polynomial _fmpz_poly_cyclotomic
#define arith_cyclotomic_polynomial fmpz_poly_cyclotomic

FLINT_DLL void _arith_cos_minpoly(fmpz * coeffs, slong d, ulong n);
FLINT_DLL void arith_cos_minpoly(fmpz_poly_t poly, ulong n);

/* Hypergeometric polynomials ************************************************/

FLINT_DLL void arith_legendre_polynomial(fmpq_poly_t poly, ulong n);

#define _arith_chebyshev_t_polynomial _fmpz_poly_chebyshev_t
#define arith_chebyshev_t_polynomial fmpz_poly_chebyshev_t
#define _arith_chebyshev_u_polynomial _fmpz_poly_chebyshev_u
#define arith_chebyshev_u_polynomial fmpz_poly_chebyshev_u

/* Swinnerton-Dyer polynomials ***********************************************/

#define _arith_swinnerton_dyer_polynomial _fmpz_poly_swinnerton_dyer
#define arith_swinnerton_dyer_polynomial fmpz_poly_swinnerton_dyer

/* Landau function ***********************************************************/

FLINT_DLL void arith_landau_function_vec(fmpz * res, slong len);

/* Dedekind sums *************************************************************/

#define arith_dedekind_sum_naive fmpq_dedekind_sum_naive
#define arith_dedekind_sum_coprime_d fmpq_dedekind_sum_coprime_d
#define arith_dedekind_sum_coprime_large fmpq_dedekind_sum_coprime_large
#define arith_dedekind_sum_coprime fmpq_dedekind_sum_coprime
#define arith_dedekind_sum fmpq_dedekind_sum

/* Exponential sums **********************************************************/

typedef struct
{
    int n;
    int prefactor;
    mp_limb_t sqrt_p;
    mp_limb_t sqrt_q;
    mp_limb_signed_t cos_p[FLINT_BITS];
    mp_limb_t cos_q[FLINT_BITS];
} trig_prod_struct;

typedef trig_prod_struct trig_prod_t[1];

ARITH_INLINE
void trig_prod_init(trig_prod_t sum)
{
    sum->n = 0;
    sum->prefactor = 1;
    sum->sqrt_p = 1;
    sum->sqrt_q = 1;
}

FLINT_DLL void arith_hrr_expsum_factored(trig_prod_t prod, mp_limb_t k, mp_limb_t n);

/* Number of partitions ******************************************************/

FLINT_DLL void arith_number_of_partitions_nmod_vec(mp_ptr res, slong len, nmod_t mod);
FLINT_DLL void arith_number_of_partitions_vec(fmpz * res, slong len);
FLINT_DLL void arith_number_of_partitions_mpfr(mpfr_t x, ulong n);
FLINT_DLL void arith_number_of_partitions(fmpz_t x, ulong n);

/* Number of sums of squares representations *********************************/

FLINT_DLL void arith_sum_of_squares(fmpz_t r, ulong k, const fmpz_t n);
FLINT_DLL void arith_sum_of_squares_vec(fmpz * r, ulong k, slong n);

#ifdef __cplusplus
}
#endif

#endif
