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

#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* MPFR extras ***************************************************************/

void mpfr_zeta_inv_euler_product(mpfr_t res, ulong s, int char_4);
void mpfr_pi_chudnovsky(mpfr_t res, mpfr_rnd_t rnd);
void mpfr_const_euler_brent_mcmillan(mpfr_t res, mpfr_rnd_t rnd);
void mpfr_zeta_ui_bsplit(mpfr_t x, ulong s, mpfr_rnd_t rnd);

/* Various arithmetic functions **********************************************/

void arith_primorial(fmpz_t res, len_t n);

void _arith_harmonic_number(fmpz_t num, fmpz_t den, len_t n);
void arith_harmonic_number(fmpq_t x, len_t n);

void arith_ramanujan_tau(fmpz_t res, const fmpz_t n);
void arith_ramanujan_tau_series(fmpz_poly_t res, len_t n);

void arith_divisors(fmpz_poly_t res, const fmpz_t n);
void arith_divisor_sigma(fmpz_t res, const fmpz_t n, ulong k);

int arith_moebius_mu(const fmpz_t n);

void arith_euler_phi(fmpz_t res, const fmpz_t n);

/* Stirling numbers **********************************************************/

void arith_stirling_number_1u(fmpz_t s, len_t n, len_t k);
void arith_stirling_number_1(fmpz_t s, len_t n, len_t k);
void arith_stirling_number_2(fmpz_t s, len_t n, len_t k);

void arith_stirling_number_1u_vec(fmpz * row, len_t n, len_t klen);
void arith_stirling_number_1_vec(fmpz * row, len_t n, len_t klen);
void arith_stirling_number_2_vec(fmpz * row, len_t n, len_t klen);

void arith_stirling_number_1u_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen);
void arith_stirling_number_1_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen);
void arith_stirling_number_2_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen);

void arith_stirling_matrix_1u(fmpz_mat_t mat);
void arith_stirling_matrix_1(fmpz_mat_t mat);
void arith_stirling_matrix_2(fmpz_mat_t mat);

/* Bell numbers **************************************************************/

#if FLINT64
#define BELL_NUMBER_TAB_SIZE 26
#else
#define BELL_NUMBER_TAB_SIZE 16
#endif

extern const mp_limb_t bell_number_tab[];

double arith_bell_number_size(ulong n);

void arith_bell_number(fmpz_t b, ulong n);
void arith_bell_number_bsplit(fmpz_t res, ulong n);
void arith_bell_number_multi_mod(fmpz_t res, ulong n);

void arith_bell_number_vec(fmpz * b, len_t n);
void arith_bell_number_vec_recursive(fmpz * b, len_t n);
void arith_bell_number_vec_multi_mod(fmpz * b, len_t n);

mp_limb_t arith_bell_number_nmod(ulong n, nmod_t mod);

void arith_bell_number_nmod_vec(mp_ptr b, len_t n, nmod_t mod);
void arith_bell_number_nmod_vec_recursive(mp_ptr b, len_t n, nmod_t mod);
void arith_bell_number_nmod_vec_series(mp_ptr b, len_t n, nmod_t mod);


/* Euler numbers *************************************************************/

#if FLINT64
#define SMALL_EULER_LIMIT 25
#else
#define SMALL_EULER_LIMIT 15
#endif

static const mp_limb_t euler_number_small[] = {
    1UL, 1UL, 5UL, 61UL, 1385UL, 50521UL, 2702765UL,
    199360981UL,
#if FLINT64
    19391512145UL, 2404879675441UL, 370371188237525UL,
    69348874393137901UL, 15514534163557086905UL
#endif
};

double arith_euler_number_size(ulong n);

void arith_euler_number_vec(fmpz * res, len_t n);

void _arith_euler_number_zeta(fmpz_t res, ulong n);
void arith_euler_number(fmpz_t res, ulong n);

void arith_euler_polynomial(fmpq_poly_t poly, ulong n);

/* Bernoulli numbers *********************************************************/

#if FLINT64
#define BERNOULLI_SMALL_NUMER_LIMIT 35
#else
#define BERNOULLI_SMALL_NUMER_LIMIT 27
#endif

static const len_t _bernoulli_numer_small[] = {
    1L, 1L, -1L, 1L, -1L, 5L, -691L, 7L, -3617L, 43867L, -174611L, 854513L,
    -236364091L, 8553103L,
#if FLINT64
    -23749461029L, 8615841276005L, -7709321041217L, 2577687858367L
#endif
};

void _arith_bernoulli_number(fmpz_t num, fmpz_t den, ulong n);
void arith_bernoulli_number(fmpq_t x, ulong n);

void _arith_bernoulli_number_vec(fmpz * num, fmpz * den, len_t n);
void arith_bernoulli_number_vec(fmpq * num, len_t n);

void arith_bernoulli_number_denom(fmpz_t den, ulong n);
double arith_bernoulli_number_size(ulong n);

void arith_bernoulli_polynomial(fmpq_poly_t poly, ulong n);

void _arith_bernoulli_number_zeta(fmpz_t num, fmpz_t den, ulong n);
void _arith_bernoulli_number_vec_multi_mod(fmpz * num, fmpz * den, len_t n);
void _arith_bernoulli_number_vec_recursive(fmpz * num, fmpz * den, len_t n);
void _arith_bernoulli_number_vec_zeta(fmpz * num, fmpz * den, len_t n);

/* Cyclotomic polynomials ****************************************************/

void _arith_cyclotomic_polynomial(fmpz * a, ulong n, mp_ptr factors,
                                        len_t num_factors, ulong phi);
void arith_cyclotomic_polynomial(fmpz_poly_t poly, ulong n);

void _arith_cos_minpoly(fmpz * coeffs, len_t d, ulong n);
void arith_cos_minpoly(fmpz_poly_t poly, ulong n);

/* Hypergeometric polynomials ************************************************/

void arith_legendre_polynomial(fmpq_poly_t poly, ulong n);
void arith_chebyshev_t_polynomial(fmpz_poly_t poly, ulong n);
void arith_chebyshev_u_polynomial(fmpz_poly_t poly, ulong n);

/* Swinnerton-Dyer polynomials ***********************************************/

void arith_swinnerton_dyer_polynomial(fmpz_poly_t poly, ulong n);

/* Landau function ***********************************************************/

void arith_landau_function_vec(fmpz * res, len_t len);

/* Dedekind sums *************************************************************/

void arith_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k);
double arith_dedekind_sum_coprime_d(double h, double k);
void arith_dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k);
void arith_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k);
void arith_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k);

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

static __inline__
void trig_prod_init(trig_prod_t sum)
{
    sum->n = 0;
    sum->prefactor = 1;
    sum->sqrt_p = 1;
    sum->sqrt_q = 1;
}

void arith_hrr_expsum_factored(trig_prod_t prod, mp_limb_t k, mp_limb_t n);

/* Number of partitions ******************************************************/

void arith_number_of_partitions_nmod_vec(mp_ptr res, len_t len, nmod_t mod);
void arith_number_of_partitions_vec(fmpz * res, len_t len);
void arith_number_of_partitions_mpfr(mpfr_t x, ulong n);
void arith_number_of_partitions(fmpz_t x, ulong n);

/* Number of sums of squares representations *********************************/

void arith_sum_of_squares(fmpz_t r, ulong k, const fmpz_t n);
void arith_sum_of_squares_vec(fmpz * r, ulong k, len_t n);

#ifdef __cplusplus
}
#endif

#endif
