/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_LLL_H
#define FMPZ_LLL_H

#ifdef FMPZ_LLL_INLINES_C
#define FMPZ_LLL_INLINE
#else
#define FMPZ_LLL_INLINE static inline
#endif

#include "d_mat.h"
#include "mpf_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

#if FLINT_BITS == 32
#define CPU_SIZE_1 31
#define FMPZ_LLL_MAX_LONG WORD(2147483647)
#else
#define CPU_SIZE_1 53
#define FMPZ_LLL_MAX_LONG WORD(9007199254740991)
#endif

#define SIZE_RED_FAILURE_THRESH 5

typedef enum
{
    GRAM,
    Z_BASIS
} rep_type;

typedef enum
{
    APPROX,
    EXACT
} gram_type;

typedef struct
{
    double delta;
    double eta;
    rep_type rt;
    gram_type gt;
} fmpz_lll_struct;

typedef fmpz_lll_struct fmpz_lll_t[1];

typedef union
{
    d_mat_t appSP;
    mpf_mat_t appSP2;
    fmpz_mat_t exactSP;
} fmpz_gram_union;

typedef fmpz_gram_union fmpz_gram_t[1];

/* Parameter manipulation  ***************************************************/

void fmpz_lll_context_init_default(fmpz_lll_t fl);

void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta,
                           rep_type rt, gram_type gt);

/* Random parameter generation  **********************************************/

void fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state);

/* The various Babai's  ******************************************************/

double fmpz_lll_heuristic_dot(const double * vec1, const double * vec2, slong len2,
       const fmpz_mat_t B, slong k, slong j, slong exp_adj);

int fmpz_lll_check_babai(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl);

int fmpz_lll_check_babai_heuristic_d(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl);

int fmpz_lll_shift(const fmpz_mat_t B);

int fmpz_lll_d(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl);

int fmpz_lll_d_heuristic(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl);

int fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U, mpf_mat_t mu, mpf_mat_t r, mpf *s,
       mpf_mat_t appB, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, mpf_t tmp, mpf_t rtmp, flint_bitcnt_t prec, const fmpz_lll_t fl);

int fmpz_lll_mpf2(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_lll_t fl);

int fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl);

int fmpz_lll_wrapper(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl);

int fmpz_lll_advance_check_babai(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl);

int fmpz_lll_advance_check_babai_heuristic_d(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl);

/* LLL with removals  ********************************************************/

int fmpz_lll_d_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_d_heuristic_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_mpf_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_wrapper_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_d_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

int fmpz_lll_wrapper_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

/* ULLL  *********************************************************************/

int fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size, const fmpz_t gs_B, const fmpz_lll_t fl);

/* LLL-reducedness ***********************************************************/

int fmpz_lll_is_reduced_d(const fmpz_mat_t B, const fmpz_lll_t fl);

int fmpz_lll_is_reduced_mpfr(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec);

int fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec);

int fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd);

int fmpz_lll_is_reduced_mpfr_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd,flint_bitcnt_t prec);

int fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec);

/* Default functions *********************************************************/

void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl);

int fmpz_lll_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl);

/* Modified ULLL  ************************************************************/

void fmpz_lll_storjohann_ulll(fmpz_mat_t FM, slong new_size, const fmpz_lll_t fl);

#ifdef __cplusplus
}
#endif

#endif

