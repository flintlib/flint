/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_MPOLY_FACTOR_H
#define FQ_ZECH_MPOLY_FACTOR_H

#ifdef FQ_ZECH_MPOLY_FACTOR_INLINES_C
#define FQ_ZECH_MPOLY_FACTOR_INLINE
#else
#define FQ_ZECH_MPOLY_FACTOR_INLINE static inline
#endif

#include "fq_zech.h"
#include "fq_zech_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fq_zech_poly_struct * coeffs;
    slong alloc;
    slong length;
} fq_zech_bpoly_struct;

typedef fq_zech_bpoly_struct fq_zech_bpoly_t[1];


typedef struct
{
    fq_zech_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_zech_tpoly_struct;

typedef fq_zech_tpoly_struct fq_zech_tpoly_t[1];


typedef struct
{
    ulong * exps;
    fq_zech_struct * coeffs;
    slong length;
    slong alloc;
} fq_zech_polyu_struct;

typedef fq_zech_polyu_struct fq_zech_polyu_t[1];


typedef struct
{
    fq_zech_poly_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
} fq_zech_polyun_struct;

typedef fq_zech_polyun_struct fq_zech_polyun_t[1];


/*****************************************************************************/

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_init(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fq_zech_bpoly_clear(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_swap(fq_zech_bpoly_t A, fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
{
    fq_zech_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_normalise(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
{
    while (A->length > 0 && fq_zech_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}

void fq_zech_bpoly_realloc(fq_zech_bpoly_t A, slong len, const fq_zech_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_fit_length(fq_zech_bpoly_t A, slong len, const fq_zech_ctx_t ctx)
{
    if (len > A->alloc)
        fq_zech_bpoly_realloc(A, len, ctx);
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_zero(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
{
    A->length = 0;
}

FQ_ZECH_MPOLY_FACTOR_INLINE
int fq_zech_bpoly_is_zero(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
{
    return A->length == 0;
}

int fq_zech_bpoly_equal(const fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_bpoly_get_coeff(fq_zech_t c, const fq_zech_bpoly_t A, slong e0, slong e1, const fq_zech_ctx_t ctx)
{
    if (e0 >= A->length)
        fq_zech_zero(c, ctx);
    else
        fq_zech_poly_get_coeff(c, A->coeffs + e0, e1, ctx);
}

FQ_ZECH_MPOLY_FACTOR_INLINE
slong fq_zech_bpoly_degree0(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
{
    return A->length - 1;
}

slong fq_zech_bpoly_degree1(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set_poly_var1(fq_zech_bpoly_t A, const fq_zech_poly_t B, const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set_poly_var0(fq_zech_bpoly_t A, const fq_zech_poly_t B, const fq_zech_ctx_t ctx);

void fq_zech_bpoly_print_pretty(const fq_zech_bpoly_t A,
                const char * var0, const char * var1, const fq_zech_ctx_t ctx);

int fq_zech_bpoly_is_canonical(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx);

int fq_zech_bpoly_fq_equal(
    const fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set_coeff_fq_zech(
    fq_zech_bpoly_t A,
    slong xi,
    slong yi,
    const fq_zech_t c,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set_fq_zech_poly_var0(
    fq_zech_bpoly_t A,
    const fq_zech_poly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set_fq_zech_poly_var1(
    fq_zech_bpoly_t A,
    const fq_zech_poly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_make_monic(
    fq_zech_bpoly_t A,
    slong order,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_mul(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_bpoly_t C,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_mul_series(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_bpoly_t C,
    slong order,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_add(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_bpoly_t C,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_one(
    fq_zech_bpoly_t A,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_sub(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_bpoly_t C,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_derivative(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_divrem_series(
    fq_zech_bpoly_t Q,
    fq_zech_bpoly_t R,
    const fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    slong order,
    const fq_zech_ctx_t ctx);

int fq_zech_bpoly_divides(
    fq_zech_bpoly_t Q,
    const fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_set(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_make_primitive(
    fq_zech_poly_t g,
    fq_zech_bpoly_t A,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_taylor_shift_var1(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    const fq_zech_t c_,
    const fq_zech_ctx_t ctx);

void fq_zech_bpoly_taylor_shift_var0(
    fq_zech_bpoly_t A,
    const fq_zech_t alpha,
    const fq_zech_ctx_t ctx);

void fq_zech_mpoly_get_fq_zech_bpoly(
    fq_zech_bpoly_t A,
    const fq_zech_mpoly_t B,
    slong varx,
    slong vary,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_fq_zech_bpoly(
    fq_zech_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_zech_bpoly_t B,
    slong varx,
    slong vary,
    const fq_zech_mpoly_ctx_t ctx);


int fq_zech_bpoly_factor_smprime(
    fq_zech_poly_t c,
    fq_zech_tpoly_t F,
    fq_zech_bpoly_t B,
    int allow_shift,
    const fq_zech_ctx_t ctx);

int fq_zech_bpoly_factor_lgprime(
    fq_zech_poly_t c,
    fq_zech_tpoly_t F,
    fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_tpoly_init(fq_zech_tpoly_t A, const fq_zech_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_tpoly_swap(fq_zech_tpoly_t A, fq_zech_tpoly_t B, const fq_zech_ctx_t ctx)
{
    fq_zech_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fq_zech_tpoly_fit_length(fq_zech_tpoly_t A, slong len, const fq_zech_ctx_t ctx);

void fq_zech_tpoly_clear(fq_zech_tpoly_t A, const fq_zech_ctx_t ctx);


/*****************************************************************************/

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyu_init(fq_zech_polyu_t A, const fq_zech_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->alloc = 0;
}

void fq_zech_polyu_clear(fq_zech_polyu_t A, const fq_zech_ctx_t ctx);

void fq_zech_polyu_realloc(fq_zech_polyu_t A, slong len, const fq_zech_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyu_fit_length(fq_zech_polyu_t A, slong len, const fq_zech_ctx_t ctx)
{
    FLINT_ASSERT(A->alloc >= 0);
    if (len > A->alloc)
        fq_zech_polyu_realloc(A, len, ctx);
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyu_swap(fq_zech_polyu_t A, fq_zech_polyu_t B, const fq_zech_ctx_t ctx)
{
    fq_zech_polyu_struct t = *B;
    *B = *A;
    *A = t;
}

void fq_zech_polyu3_print_pretty(
    const fq_zech_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fq_zech_ctx_t ctx);

void fq_zech_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const fq_zech_polyu_t A);

int fq_zech_polyu_is_canonical(
    const fq_zech_polyu_t A,
    const fq_zech_ctx_t ctx);

/*****************************************************************************/

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyun_init(fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->alloc = 0;
}

void fq_zech_polyun_clear(fq_zech_polyun_t A, const fq_zech_ctx_t ctx);

void fq_zech_polyun_realloc(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyun_fit_length(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx)
{
    if (len > A->alloc)
        fq_zech_polyun_realloc(A, len, ctx);
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_polyun_swap(fq_zech_polyun_t A, fq_zech_polyun_t B, const fq_zech_ctx_t ctx)
{
    fq_zech_polyun_struct t = *B;
    *B = *A;
    *A = t;
}

void fq_zech_polyu2n_print_pretty(
    const fq_zech_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_zech_ctx_t ctx);

void fq_zech_polyu3n_print_pretty(
    const fq_zech_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fq_zech_ctx_t ctx);

int fq_zech_polyun_is_canonical(
    const fq_zech_polyun_t A,
    const fq_zech_ctx_t ctx);

/*****************************************************************************/

int fq_zech_mpoly_is_fq_zech_poly(
    const fq_zech_mpoly_t A,
    slong var,
    const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_get_fq_zech_poly(
    fq_zech_poly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_set_fq_zech_poly(
    fq_zech_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_zech_struct * Bcoeffs,
    slong Blen,
    slong var,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_fq_zech_poly(
    fq_zech_mpoly_t A,
    const fq_zech_poly_t B,
    slong var,
    const fq_zech_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct {
    fq_zech_t constant;
    fq_zech_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fq_zech_mpoly_factor_struct;

typedef fq_zech_mpoly_factor_struct fq_zech_mpoly_factor_t[1];

void fq_zech_mpoly_factor_init(fq_zech_mpoly_factor_t f,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_realloc(fq_zech_mpoly_factor_t f,
                                   slong alloc, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_fit_length(fq_zech_mpoly_factor_t f,
                                     slong len, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_clear(fq_zech_mpoly_factor_t f,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_set(fq_zech_mpoly_factor_t a,
                const fq_zech_mpoly_factor_t b, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_print_pretty(const fq_zech_mpoly_factor_t f,
                            const char ** vars, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_append_ui(fq_zech_mpoly_factor_t f,
              const fq_zech_mpoly_t A, ulong e, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_factor_append_fmpz(fq_zech_mpoly_factor_t f,
       const fq_zech_mpoly_t A, const fmpz_t e, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_factor_squarefree(fq_zech_mpoly_factor_t f,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_factor(fq_zech_mpoly_factor_t f,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_mpoly_factor_swap(fq_zech_mpoly_factor_t A,
                       fq_zech_mpoly_factor_t B, const fq_zech_mpoly_ctx_t ctx)
{
   fq_zech_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_mpoly_factor_one(fq_zech_mpoly_factor_t a,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
	fq_zech_one(a->constant, ctx->fqctx);
	a->num = 0;
}

int fq_zech_mpoly_factor_expand(fq_zech_mpoly_t A,
                const fq_zech_mpoly_factor_t f, const fq_zech_mpoly_ctx_t ctx);


FQ_ZECH_MPOLY_FACTOR_INLINE
int fq_zech_mpoly_factor_matches(const fq_zech_mpoly_t a, const fq_zech_mpoly_factor_t f, const fq_zech_mpoly_ctx_t ctx)
{
    int matches;
    fq_zech_mpoly_t t;
    fq_zech_mpoly_init(t, ctx);
    fq_zech_mpoly_factor_expand(t, f, ctx);
    matches = fq_zech_mpoly_equal(t, a, ctx);
    fq_zech_mpoly_clear(t, ctx);
    return matches;
}

void _fq_zech_mpoly_get_lead0(
    fq_zech_mpoly_t c,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_set_lead0(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_t c,
    const fq_zech_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fq_zech_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_zech_mpolyv_struct;

typedef fq_zech_mpolyv_struct fq_zech_mpolyv_t[1];

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_mpolyv_init(fq_zech_mpolyv_t A, const fq_zech_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_ZECH_MPOLY_FACTOR_INLINE
void fq_zech_mpolyv_swap(fq_zech_mpolyv_t A, fq_zech_mpolyv_t B,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
   fq_zech_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_zech_mpolyv_clear(fq_zech_mpolyv_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyv_print_pretty(const fq_zech_mpolyv_t poly,
                               const char ** x, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyv_fit_length(fq_zech_mpolyv_t A, slong length,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyv_set_coeff(
    fq_zech_mpolyv_t A,
    slong i,
    fq_zech_mpoly_t c,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_to_mpolyv(
    fq_zech_mpolyv_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_t xalpha,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_from_mpolyv(
    fq_zech_mpoly_t A,
    const fq_zech_mpolyv_t B,
    const fq_zech_mpoly_t xalpha,
    const fq_zech_mpoly_ctx_t ctx);

/*****************************************************************************/

int fq_zech_mpoly_univar_content_mpoly(
    fq_zech_mpoly_t g,
    const fq_zech_mpoly_univar_t A,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_univar_divexact_mpoly(
    fq_zech_mpoly_univar_t A,
    const fq_zech_mpoly_t b,
    const fq_zech_mpoly_ctx_t ctx);

/*****************************************************************************/

int fq_zech_mpoly_factor_lcc_wang(
    fq_zech_mpoly_struct * lc_divs,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_poly_t Auc,
    const fq_zech_bpoly_struct * Auf,
    slong r,
    const fq_zech_poly_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_factor_irred_smprime_zassenhaus(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_zech_mpoly_factor_irred_lgprime_zassenhaus(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_zech_mpoly_factor_irred_smprime_wang(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_zech_mpoly_factor_irred_lgprime_wang(
    fq_zech_mpolyv_t Af,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_zech_mpoly_factor_irred_smprime_zippel(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_zech_mpoly_factor_irred_lgprime_zippel(
    fq_zech_mpolyv_t Af,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    fq_zech_poly_struct * inv_prod_dbetas;
    fq_zech_mpoly_struct * inv_prod_dbetas_mvar;
    fq_zech_poly_struct * dbetas;
    fq_zech_mpoly_struct * dbetas_mvar;
    fq_zech_mpoly_struct * prod_mbetas;
    fq_zech_mpolyv_struct * prod_mbetas_coeffs;
    fq_zech_mpoly_struct * mbetas;
    fq_zech_mpoly_struct * deltas;
    fq_zech_mpoly_struct * xalpha;
    fq_zech_mpoly_struct * q;
    fq_zech_mpoly_struct * qt;
    fq_zech_mpoly_struct * newt;
    fq_zech_mpolyv_struct * delta_coeffs;
    fq_zech_mpoly_t T;
    fq_zech_mpoly_t Q;
    fq_zech_mpoly_t R;
} fq_zech_mpoly_pfrac_struct;

typedef fq_zech_mpoly_pfrac_struct fq_zech_mpoly_pfrac_t[1];


int fq_zech_mpoly_pfrac_init(
    fq_zech_mpoly_pfrac_t I,
    flint_bitcnt_t bits,
    slong l, slong r,
    const fq_zech_mpoly_struct * betas,
    const fq_zech_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_pfrac_clear(
    fq_zech_mpoly_pfrac_t I,
    const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_pfrac(
    slong r,
    fq_zech_mpoly_t t,
    const slong * deg,
    fq_zech_mpoly_pfrac_t I,
    const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_hlift(
    slong m,
    fq_zech_mpoly_struct * f, /* length r */
    slong r,
    const fq_zech_struct * alpha,
    const fq_zech_mpoly_t A,
    const slong * degs,
    const fq_zech_mpoly_ctx_t ctx);

int fq_zech_bpoly_hlift2(
    fq_zech_bpoly_t A, /* clobbered (shifted by alpha) */
    fq_zech_bpoly_t B0,
    fq_zech_bpoly_t B1,
    const fq_zech_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_zech_ctx_t ctx);

int fq_zech_bpoly_hlift(
    slong r,
    fq_zech_bpoly_t A, /* clobbered (shifted by alpha) */
    fq_zech_bpoly_struct * B,
    const fq_zech_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_zech_ctx_t ctx);

int fq_zech_polyu3_hlift(
    slong r,
    fq_zech_polyun_struct * BB,
    fq_zech_polyu_t A,
    fq_zech_polyu_struct * B,
    const fq_zech_t beta,
    slong degree_inner, /* required degree in x */
    const fq_zech_ctx_t ctx);

int fq_zech_mpoly_factor_algo(fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx, unsigned int algo);

int fq_zech_mpoly_factor_zassenhaus(fq_zech_mpoly_factor_t f,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_factor_wang(fq_zech_mpoly_factor_t f,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_factor_zippel(fq_zech_mpoly_factor_t f,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_poly_product_roots_fq_zech(
    fq_zech_poly_t master,
    const fq_zech_struct * monomials,
    slong mlength,
    const fq_zech_ctx_t ctx);

void _fq_zech_mpoly_monomial_evals(
    fq_zech_struct * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_zech_struct * alpha,
    slong vstart,
    const fq_zech_mpoly_ctx_t ctx);

int _fq_zech_mpoly_eval_rest_fq_zech_poly(
    fq_zech_poly_struct * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const fq_zech_struct * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const fq_zech_poly_struct * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    const fq_zech_ctx_t ctx);

void _fq_zech_mpoly_eval_to_bpoly(
    fq_zech_bpoly_t E,
    const fq_zech_mpoly_t A,
    const fq_zech_poly_struct * alphabetas,
    const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(
    fq_zech_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_zech_bpoly_t B,
    slong var,
    const fq_zech_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

