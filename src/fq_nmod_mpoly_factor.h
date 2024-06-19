/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_FACTOR_H
#define FQ_NMOD_MPOLY_FACTOR_H

#ifdef FQ_NMOD_MPOLY_FACTOR_INLINES_C
#define FQ_NMOD_MPOLY_FACTOR_INLINE
#else
#define FQ_NMOD_MPOLY_FACTOR_INLINE static inline
#endif

#include "fq_nmod_mpoly.h"

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

void fq_nmod_mpoly_factor_init(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_realloc(fq_nmod_mpoly_factor_t f,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_fit_length(fq_nmod_mpoly_factor_t f,
                                     slong len, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_clear(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
slong fq_nmod_mpoly_factor_length(const fq_nmod_mpoly_factor_t f,
                                                 const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    return f->num;
}

void fq_nmod_mpoly_factor_get_constant_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx);

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

slong fq_nmod_mpoly_factor_get_exp_si(fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx));

void fq_nmod_mpoly_factor_set(fq_nmod_mpoly_factor_t a,
                const fq_nmod_mpoly_factor_t b, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_print_pretty(const fq_nmod_mpoly_factor_t f,
                            const char ** vars, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_append_ui(fq_nmod_mpoly_factor_t f,
              const fq_nmod_mpoly_t A, ulong e, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_append_fmpz(fq_nmod_mpoly_factor_t f,
       const fq_nmod_mpoly_t A, const fmpz_t e, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_content(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_squarefree(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_separable(fq_nmod_mpoly_factor_t f,
              const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, int sep);

int fq_nmod_mpoly_factor(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_factor_sort(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_cmp(const fq_nmod_mpoly_factor_t A,
                const fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_swap(fq_nmod_mpoly_factor_t A,
                       fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
   fq_nmod_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpoly_factor_one(fq_nmod_mpoly_factor_t a, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_expand(fq_nmod_mpoly_t A,
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

void _fq_nmod_mpoly_get_lead0(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_set_lead0(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

void n_fq_bpoly_mul(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_mul_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_add(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_sub(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_divrem_series(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx);

int n_fq_bpoly_divides(
    n_bpoly_t Q,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_make_primitive(
    n_poly_t g,
    n_bpoly_t A,
    const fq_nmod_ctx_t ctx);

void fq_nmod_mpoly_get_n_fq_bpoly(
    n_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_n_fq_bpoly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);


int n_fq_bpoly_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    const fq_nmod_ctx_t ctx);

int n_fq_bpoly_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

void n_polyu3_fq_print_pretty(
    const n_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fq_nmod_ctx_t ctx);

int n_polyu_fq_is_canonical(
    const n_polyu_t A,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

void n_polyu2n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_nmod_ctx_t ctx);

void n_polyu3n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fq_nmod_ctx_t ctx);

int n_polyun_fq_is_canonical(
    const n_polyun_t A,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fq_nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_mpolyv_struct;

typedef fq_nmod_mpolyv_struct fq_nmod_mpolyv_t[1];

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_init(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_swap(fq_nmod_mpolyv_t A, fq_nmod_mpolyv_t B,
                                                 const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
   fq_nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyv_clear(fq_nmod_mpolyv_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyv_print_pretty(const fq_nmod_mpolyv_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyv_fit_length(fq_nmod_mpolyv_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyv_set_coeff(fq_nmod_mpolyv_t A, slong i,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_mpolyv(fq_nmod_mpolyv_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t xalpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_from_mpolyv(fq_nmod_mpoly_t A,
                 flint_bitcnt_t Abits, const fq_nmod_mpolyv_t B,
                 const fq_nmod_mpoly_t xalpha, const fq_nmod_mpoly_ctx_t ctx);

int _fq_nmod_mpoly_vec_content_mpoly(fq_nmod_mpoly_t g,
                                  const fq_nmod_mpoly_struct * A, slong Alen,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_vec_divexact_mpoly(fq_nmod_mpoly_struct * A,
           slong Alen, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_vec_mul_mpoly(fq_nmod_mpoly_struct * A,
           slong Alen, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

int _fq_nmod_mpoly_factor_separable(fq_nmod_mpoly_factor_t f,
              const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, int sep);

int fq_nmod_mpoly_factor_lcc_wang(
    fq_nmod_mpoly_struct * lc_divs,
    const fq_nmod_mpoly_factor_t lcAfac,
    const n_poly_t Auc,
    const n_bpoly_struct * Auf,
    slong r,
    const n_poly_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_irred_smprime_zassenhaus(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_nmod_mpoly_factor_irred_lgprime_zassenhaus(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_nmod_mpoly_factor_irred_smprime_wang(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_nmod_mpoly_factor_irred_lgprime_wang(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_nmod_mpoly_factor_irred_smprime_zippel(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fq_nmod_mpoly_factor_irred_lgprime_zippel(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

void fq_nmod_mpoly_compression_do(fq_nmod_mpoly_t L,
               const fq_nmod_mpoly_ctx_t Lctx, ulong * Acoeffs, slong Alen,
                                                        mpoly_compression_t M);

void fq_nmod_mpoly_compression_undo(fq_nmod_mpoly_t A,
       flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t Actx, fq_nmod_mpoly_t L,
                        const fq_nmod_mpoly_ctx_t Lctx, mpoly_compression_t M);

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
    fq_nmod_mpoly_geobucket_struct * G;
    fq_nmod_mpoly_struct * qt;
    fq_nmod_mpoly_struct * newt;
    fq_nmod_mpolyv_struct * delta_coeffs;
    fq_nmod_mpoly_t T;
    fq_nmod_mpoly_t Q;
    fq_nmod_mpoly_t R;
} fq_nmod_mpoly_pfrac_struct;

typedef fq_nmod_mpoly_pfrac_struct fq_nmod_mpoly_pfrac_t[1];


int fq_nmod_mpoly_pfrac_init(
    fq_nmod_mpoly_pfrac_t Iv,
    flint_bitcnt_t bits,
    slong l, slong r,
    const fq_nmod_mpoly_struct * betas,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_pfrac_clear(
    fq_nmod_mpoly_pfrac_t Iv,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_pfrac(
    slong r,
    fq_nmod_mpoly_t t,
    const slong * deg,
    fq_nmod_mpoly_pfrac_t Iv,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_hlift(
    slong m,
    fq_nmod_mpoly_struct * f, /* length r */
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx);

int n_fq_bpoly_hlift2_cubic(
    n_fq_bpoly_t A, /* clobbered (shifted by alpha) */
    n_fq_bpoly_t B0,
    n_fq_bpoly_t B1,
    const fq_nmod_t alpha_,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    nmod_eval_interp_t E,
    n_poly_bpoly_stack_t St);

int n_fq_bpoly_hlift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t St);

int n_fq_bpoly_hlift_cubic(
    n_fq_bpoly_t A, /* clobbered (shifted by alpha) */
    n_fq_bpoly_t B0,
    n_fq_bpoly_t B1,
    const fq_nmod_t alpha_,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    nmod_eval_interp_t E,
    n_poly_bpoly_stack_t St);

int n_fq_bpoly_hlift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t FLINT_UNUSED(St));

int n_fq_polyu3_hlift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    const fq_nmod_t beta,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t St);

int fq_nmod_mpoly_factor_algo(fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, unsigned int algo);

int fq_nmod_mpoly_factor_wang(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_zassenhaus(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_factor_zippel(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int _fq_nmod_mpoly_eval_rest_n_fq_poly(n_poly_struct * E,
    slong * starts, slong * ends, slong * stops, ulong * es,
    const ulong * Acoeffs, const ulong * Aexps, slong Alen, slong var,
    const n_fq_poly_struct * alphas, const slong * offsets, const slong * shifts,
                    slong N, ulong mask, slong nvars, const fq_nmod_ctx_t ctx);

void _fq_nmod_mpoly_eval_rest_to_n_fq_bpoly(n_bpoly_t E,
                    const fq_nmod_mpoly_t A, const n_poly_struct * alphabetas,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_set_n_fq_bpoly_gen1_zero(fq_nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const n_bpoly_t B, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
