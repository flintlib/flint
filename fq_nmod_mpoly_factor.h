/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_FACTOR_H
#define FQ_NMOD_MPOLY_FACTOR_H

#ifdef FQ_NMOD_MPOLY_FACTOR_INLINES_C
#define FQ_NMOD_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "fq_nmod_mpoly.h"
#include "nmod_mpoly_factor.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*****************************************************************************/

typedef struct {
    fq_nmod_t constant;
    fq_nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fq_nmod_mpoly_factor_struct;

typedef fq_nmod_mpoly_factor_struct fq_nmod_mpoly_factor_t[1];

FLINT_DLL void fq_nmod_mpoly_factor_init(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_realloc(fq_nmod_mpoly_factor_t f,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_fit_length(fq_nmod_mpoly_factor_t f,
                                     slong len, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_clear(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
slong fq_nmod_mpoly_factor_length(const fq_nmod_mpoly_factor_t f,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return f->num;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_get_constant_fq_nmod(fq_nmod_t c,
                 const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_set(c, f->constant, ctx->fqctx);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_get_base(fq_nmod_mpoly_t p,
        const fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fq_nmod_mpoly_set(p, f->poly + i, ctx);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_swap_base(fq_nmod_mpoly_t p,
        const fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fq_nmod_mpoly_swap(p, f->poly + i, ctx);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
slong fq_nmod_mpoly_factor_get_exp_si(fq_nmod_mpoly_factor_t f,
                                        slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}


FLINT_DLL void fq_nmod_mpoly_factor_set(fq_nmod_mpoly_factor_t a,
                const fq_nmod_mpoly_factor_t b, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_print_pretty(const fq_nmod_mpoly_factor_t f,
                            const char ** vars, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_append_ui(fq_nmod_mpoly_factor_t f,
              const fq_nmod_mpoly_t A, ulong e, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_append_fmpz(fq_nmod_mpoly_factor_t f,
       const fq_nmod_mpoly_t A, const fmpz_t e, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_squarefree(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_sort(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_cmp(const fq_nmod_mpoly_factor_t A,
                const fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_swap(fq_nmod_mpoly_factor_t A,
                       fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_one(fq_nmod_mpoly_factor_t a,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
	fq_nmod_one(a->constant, ctx->fqctx);
	a->num = 0;
}

FLINT_DLL int fq_nmod_mpoly_factor_expand(fq_nmod_mpoly_t A,
                const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx);


FQ_NMOD_MPOLY_FACTOR_INLINE
int fq_nmod_mpoly_factor_matches(const fq_nmod_mpoly_t a,
                 const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
{
    int matches;
    fq_nmod_mpoly_t t;
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_factor_expand(t, f, ctx);
    matches = fq_nmod_mpoly_equal(t, a, ctx);
    fq_nmod_mpoly_clear(t, ctx);
    return matches;
}

FLINT_DLL void _fq_nmod_mpoly_get_lead0(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_lead0(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL void bad_n_fq_embed_sm_to_lg(
    mp_limb_t * out,            /* element of lgctx */
    const n_poly_t in,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb);

FLINT_DLL void bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(
    fq_nmod_t out,            /* element of lgctx */
    const n_poly_t in,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb);

FLINT_DLL void bad_n_fq_embed_lg_to_sm(
    n_poly_t out,  /* poly over smctx */
    const mp_limb_t * in,  /* element of lgctx */
    const bad_fq_nmod_embed_t emb);

/*****************************************************************************/

FLINT_DLL void n_bpoly_fq_print_pretty(const n_bpoly_t A,
                const char * var0, const char * var1, const fq_nmod_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_is_canonical(const n_bpoly_t A, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_one(n_bpoly_t A, const fq_nmod_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_equal(
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_coeff_n_fq(
    n_bpoly_t A,
    slong xi,
    slong yi,
    const mp_limb_t *c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_coeff_fq_nmod(
    n_bpoly_t A,
    slong xi,
    slong yi,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_fq_nmod_poly_var0(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_fq_nmod_poly_var1(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_make_monic(
    n_bpoly_t A,
    slong order,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_mul(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_mul_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_add(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_one(n_bpoly_t A, const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_sub(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_derivative(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_divrem_series(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_divides(
    n_bpoly_t Q,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_make_primitive(
    n_poly_t g,
    n_bpoly_t A,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_taylor_shift_var1_fq_nmod(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_t c_,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_taylor_shift_var0_fq_nmod(
    n_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_n_bpoly_fq(
    n_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_n_bpoly_fq(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);


FLINT_DLL int n_bpoly_fq_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

FLINT_DLL void n_polyu3_fq_print_pretty(
    const n_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_polyu_fq_is_canonical(
    const n_polyu_t A,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL void n_polyu2n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_polyu3n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_polyun_fq_is_canonical(
    const n_polyun_t A,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int fq_nmod_mpoly_is_fq_nmod_poly(
    const fq_nmod_mpoly_t A,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_get_fq_nmod_poly(
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_fit_length_set_bits(
    fq_nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_struct * Bcoeffs,
    slong Blen,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    const fq_nmod_poly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fq_nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_mpolyv_struct;

typedef fq_nmod_mpolyv_struct fq_nmod_mpolyv_t[1];

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_init(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_swap(fq_nmod_mpolyv_t A, fq_nmod_mpolyv_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void fq_nmod_mpolyv_clear(fq_nmod_mpolyv_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyv_print_pretty(const fq_nmod_mpolyv_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyv_fit_length(fq_nmod_mpolyv_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyv_set_coeff(
    fq_nmod_mpolyv_t A,
    slong i,
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_to_mpolyv(
    fq_nmod_mpolyv_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t xalpha,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_from_mpolyv(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpolyv_t B,
    const fq_nmod_mpoly_t xalpha,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int _fq_nmod_mpoly_vec_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_struct * A,
    slong Alen,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

int fq_nmod_mpoly_factor_lcc_wang(
    fq_nmod_mpoly_struct * lc_divs,
    const fq_nmod_mpoly_factor_t lcAfac,
    const n_poly_t Auc,
    const n_bpoly_struct * Auf,
    slong r,
    const n_poly_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_irred_smprime_zassenhaus(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_lgprime_zassenhaus(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_smprime_wang(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_lgprime_wang(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_smprime_zippel(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_lgprime_zippel(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    fq_nmod_poly_struct * inv_prod_dbetas;
    fq_nmod_mpoly_struct * inv_prod_dbetas_mvar;
    fq_nmod_poly_struct * dbetas;
    fq_nmod_mpoly_struct * dbetas_mvar;
    fq_nmod_mpoly_struct * prod_mbetas;
    fq_nmod_mpolyv_struct * prod_mbetas_coeffs;
    fq_nmod_mpoly_struct * mbetas;
    fq_nmod_mpoly_struct * deltas;
    fq_nmod_mpoly_struct * xalpha;
    fq_nmod_mpoly_struct * q;
    fq_nmod_mpoly_struct * qt;
    fq_nmod_mpoly_struct * newt;
    fq_nmod_mpolyv_struct * delta_coeffs;
    fq_nmod_mpoly_t T;
    fq_nmod_mpoly_t Q;
    fq_nmod_mpoly_t R;
} fq_nmod_mpoly_pfrac_struct;

typedef fq_nmod_mpoly_pfrac_struct fq_nmod_mpoly_pfrac_t[1];


FLINT_DLL int fq_nmod_mpoly_pfrac_init(
    fq_nmod_mpoly_pfrac_t I,
    flint_bitcnt_t bits,
    slong l, slong r,
    const fq_nmod_mpoly_struct * betas,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_pfrac_clear(
    fq_nmod_mpoly_pfrac_t I,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_pfrac(
    slong r,
    fq_nmod_mpoly_t t,
    const slong * deg,
    fq_nmod_mpoly_pfrac_t I,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_hlift(
    slong m,
    fq_nmod_mpoly_struct * f, /* length r */
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_hlift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_bpoly_fq_hlift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_polyu3_fq_hlift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    const fq_nmod_t beta,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_product_roots_n_fq(
    n_poly_t master,
    const mp_limb_t * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_monomial_evals(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * alpha,
    slong vstart,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_wang(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_zassenhaus(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_zippel(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int _fq_nmod_mpoly_eval_rest_n_poly_fq(n_poly_struct * E,
    slong * starts, slong * ends, slong * stops, ulong * es,
    const fq_nmod_struct * Acoeffs, const ulong * Aexps, slong Alen, slong var,
    const n_poly_struct * alphas, const slong * offsets, const slong * shifts, 
                    slong N, ulong mask, slong nvars, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_eval_rest_to_n_bpoly_fq(n_bpoly_t E,
                    const fq_nmod_mpoly_t A, const n_poly_struct * alphabetas,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_n_bpoly_fq_var1_zero(fq_nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const n_bpoly_t B, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

