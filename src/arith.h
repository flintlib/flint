/*
    Copyright (C) 2010-2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARITH_H
#define ARITH_H

#include "fmpq_types.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Various arithmetic functions **********************************************/

void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n);
void arith_harmonic_number(fmpq_t x, slong n);

void arith_ramanujan_tau(fmpz_t res, const fmpz_t n);
void arith_ramanujan_tau_series(fmpz_poly_t res, slong n);

void arith_divisors(fmpz_poly_t res, const fmpz_t n);

/* Stirling numbers **********************************************************/

void arith_stirling_number_1u(fmpz_t s, ulong n, ulong k);
void arith_stirling_number_1(fmpz_t s, ulong n, ulong k);
void arith_stirling_number_2(fmpz_t s, ulong n, ulong k);

void arith_stirling_number_1u_vec(fmpz * row, ulong n, slong klen);
void arith_stirling_number_1_vec(fmpz * row, ulong n, slong klen);
void arith_stirling_number_2_vec(fmpz * row, ulong n, slong klen);

void arith_stirling_number_1u_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);
void arith_stirling_number_1_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);
void arith_stirling_number_2_vec_next(fmpz * row,
        const fmpz * prev, slong n, slong klen);

void arith_stirling_matrix_1u(fmpz_mat_t mat);
void arith_stirling_matrix_1(fmpz_mat_t mat);
void arith_stirling_matrix_2(fmpz_mat_t mat);

/* Bell numbers **************************************************************/

#if FLINT64
#define BELL_NUMBER_TAB_SIZE 26
#else
#define BELL_NUMBER_TAB_SIZE 16
#endif

FLINT_DLL extern const ulong bell_number_tab[];

double arith_bell_number_size(ulong n);

void arith_bell_number(fmpz_t b, ulong n);
void arith_bell_number_dobinski(fmpz_t res, ulong n);
void arith_bell_number_multi_mod(fmpz_t res, ulong n);

void arith_bell_number_vec(fmpz * b, slong n);
void arith_bell_number_vec_recursive(fmpz * b, slong n);
void arith_bell_number_vec_multi_mod(fmpz * b, slong n);

ulong arith_bell_number_nmod(ulong n, nmod_t mod);

void arith_bell_number_nmod_vec(nn_ptr b, slong n, nmod_t mod);
void arith_bell_number_nmod_vec_recursive(nn_ptr b, slong n, nmod_t mod);
int arith_bell_number_nmod_vec_series(nn_ptr b, slong n, nmod_t mod);
void arith_bell_number_nmod_vec_ogf(nn_ptr res, slong len, nmod_t mod);


/* Euler numbers *************************************************************/

#if FLINT64
#define SMALL_EULER_LIMIT 25
#else
#define SMALL_EULER_LIMIT 15
#endif

static const ulong euler_number_small[] = {
    UWORD(1), UWORD(1), UWORD(5), UWORD(61), UWORD(1385), UWORD(50521), UWORD(2702765),
    UWORD(199360981),
#if FLINT64
    UWORD(19391512145), UWORD(2404879675441), UWORD(370371188237525),
    UWORD(69348874393137901), UWORD(15514534163557086905)
#endif
};

double arith_euler_number_size(ulong n);
void arith_euler_number_vec(fmpz * res, slong n);
void arith_euler_number(fmpz_t res, ulong n);

void arith_euler_polynomial(fmpq_poly_t poly, ulong n);

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

void _arith_bernoulli_number(fmpz_t num, fmpz_t den, ulong n);
void arith_bernoulli_number(fmpq_t x, ulong n);

void _arith_bernoulli_number_vec(fmpz * num, fmpz * den, slong n);
void arith_bernoulli_number_vec(fmpq * num, slong n);

void arith_bernoulli_number_denom(fmpz_t den, ulong n);
double arith_bernoulli_number_size(ulong n);

void arith_bernoulli_polynomial(fmpq_poly_t poly, ulong n);

void _arith_bernoulli_number_vec_multi_mod(fmpz * num, fmpz * den, slong n);
void _arith_bernoulli_number_vec_recursive(fmpz * num, fmpz * den, slong n);

/* Landau function ***********************************************************/

void arith_landau_function_vec(fmpz * res, slong len);

/* Exponential sums **********************************************************/

typedef struct
{
    int n;
    int prefactor;
    ulong sqrt_p;
    ulong sqrt_q;
    slong cos_p[FLINT_BITS];
    ulong cos_q[FLINT_BITS];
} trig_prod_struct;

typedef trig_prod_struct trig_prod_t[1];

FLINT_FORCE_INLINE
void trig_prod_init(trig_prod_t sum)
{
    sum->n = 0;
    sum->prefactor = 1;
    sum->sqrt_p = 1;
    sum->sqrt_q = 1;
}

void arith_hrr_expsum_factored(trig_prod_t prod, ulong k, ulong n);

/* Number of partitions ******************************************************/

FLINT_DLL extern const unsigned int partitions_lookup[128];

void arith_number_of_partitions_nmod_vec(nn_ptr res, slong len, nmod_t mod);
void arith_number_of_partitions_vec(fmpz * res, slong len);
void arith_number_of_partitions(fmpz_t x, ulong n);

/* Number of sums of squares representations *********************************/

void arith_sum_of_squares(fmpz_t r, ulong k, const fmpz_t n);
void arith_sum_of_squares_vec(fmpz * r, ulong k, slong n);

/* to be deprecated **********************************************************/

#define arith_primorial fmpz_primorial
#define arith_divisor_sigma fmpz_divisor_sigma
#define arith_moebius_mu fmpz_moebius_mu
#define arith_euler_phi fmpz_euler_phi
#define _arith_cyclotomic_polynomial _fmpz_poly_cyclotomic
#define arith_cyclotomic_polynomial fmpz_poly_cyclotomic
#define _arith_chebyshev_t_polynomial _fmpz_poly_chebyshev_t
#define arith_chebyshev_t_polynomial fmpz_poly_chebyshev_t
#define _arith_chebyshev_u_polynomial _fmpz_poly_chebyshev_u
#define arith_chebyshev_u_polynomial fmpz_poly_chebyshev_u
#define arith_legendre_polynomial fmpq_poly_legendre_p
#define _arith_swinnerton_dyer_polynomial _fmpz_poly_swinnerton_dyer
#define arith_swinnerton_dyer_polynomial fmpz_poly_swinnerton_dyer
#define arith_dedekind_sum_naive fmpq_dedekind_sum_naive
#define arith_dedekind_sum fmpq_dedekind_sum

#ifdef __cplusplus
}
#endif

#endif
