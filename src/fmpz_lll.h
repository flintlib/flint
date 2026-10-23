/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_LLL_H
#define FMPZ_LLL_H

#ifdef FMPZ_LLL_INLINES_C
#define FMPZ_LLL_INLINE
#else
#define FMPZ_LLL_INLINE static inline
#endif

#include "fmpz_types.h"
#include "d_mat.h"
#include "gr_types.h"

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
    gr_mat_t appSP2;
    fmpz_mat_t exactSP;
} fmpz_gram_union;

typedef fmpz_gram_union fmpz_gram_t[1];

/* Packed basis matrix  ******************************************************/

typedef struct
{
    nn_ptr entries;     /* d * n * m limbs */
    nn_ptr * rows;      /* row pointers, permuted by row moves */
    slong * bits;       /* upper bound on the bit length of the entries of each row */
    slong d;            /* rows */
    slong n;            /* columns */
    slong m;            /* limbs per entry (two's complement) */
}
fmpz_lll_packed_struct;

typedef fmpz_lll_packed_struct fmpz_lll_packed_t[1];

/*
    Only pack if the entries fit in this many limbs. In the multiprecision
    LLL the Gram matrix is computed exactly from the packed rows, which is
    only competitive with approximate dot products for small entries.
*/
#define FMPZ_LLL_PACKED_MAX_LIMBS 32
#define FMPZ_LLL_PACKED_MAX_LIMBS_MPF FLINT_MPN_DOT_TAB_N

/* do not pack if the packed matrix would need more limbs than this */
#define FMPZ_LLL_PACKED_MAX_SIZE (WORD(1) << 26)

slong fmpz_lll_packed_limbs(const fmpz_mat_t B);
void fmpz_lll_packed_init(fmpz_lll_packed_t P, slong d, slong n, slong m);
void fmpz_lll_packed_clear(fmpz_lll_packed_t P);
void fmpz_lll_packed_set_fmpz_mat(fmpz_lll_packed_t P, const fmpz_mat_t B);
void fmpz_lll_packed_get_fmpz_mat(fmpz_mat_t B, const fmpz_lll_packed_t P);
void fmpz_lll_packed_tighten(fmpz_lll_packed_t P, slong i);
void fmpz_lll_packed_grow(fmpz_lll_packed_t P);
void fmpz_lll_packed_maybe_shrink(fmpz_lll_packed_t P);
void fmpz_lll_packed_move_row(fmpz_lll_packed_t P, slong i, slong j);
void fmpz_lll_packed_row_sub(fmpz_lll_packed_t P, slong i, slong j);
void fmpz_lll_packed_row_add(fmpz_lll_packed_t P, slong i, slong j);
void fmpz_lll_packed_row_submul_si(fmpz_lll_packed_t P, slong i, slong j, slong x);
void fmpz_lll_packed_row_submul_fmpz(fmpz_lll_packed_t P, slong i, slong j, const fmpz_t x);
void fmpz_lll_packed_row_submul_si_2exp(fmpz_lll_packed_t P, slong i, slong j, slong x, ulong e);
slong fmpz_lll_packed_get_d_vec_2exp(double * appv, fmpz_lll_packed_t P, slong i);
void fmpz_lll_packed_dot(fmpz_t res, const fmpz_lll_packed_t P, slong i, slong j, slong len);
double fmpz_lll_packed_heuristic_dot(const double * vec1, const double * vec2, slong len2,
       const fmpz_lll_packed_t P, slong k, slong j, slong exp_adj);

int fmpz_lll_check_babai_packed(int kappa, fmpz_lll_packed_t P, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl, int heuristic);
int fmpz_lll_advance_check_babai_packed(int cur_kappa, int kappa, fmpz_lll_packed_t P, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s,
       d_mat_t appB, int *expo, fmpz_gram_t A,
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl, int heuristic);
int fmpz_lll_check_babai_heuristic_packed(int kappa, fmpz_lll_packed_t P, fmpz_mat_t U,
       gr_mat_t mu, gr_mat_t r, gr_ptr s, fmpz_gram_t A, int a, int zeros,
       int kappamax, int n, gr_ptr tmp, gr_ptr rtmp, gr_ctx_t ctx, const fmpz_lll_t fl);

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

int fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U,
                               gr_mat_t mu, gr_mat_t r, gr_ptr s,
                               gr_mat_t appB, fmpz_gram_t A, int a, int zeros,
                               int kappamax, int n, gr_ptr tmp, gr_ptr rtmp,
                               gr_ctx_t ctx, const fmpz_lll_t fl);

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
