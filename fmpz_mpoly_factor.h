/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_FACTOR_H
#define FMPZ_MPOLY_FACTOR_H

#ifdef FMPZ_MPOLY_FACTOR_INLINES_C
#define FMPZ_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "fmpz_mpoly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_factor.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*****************************************************************************/

FLINT_DLL void fmpz_mpoly_fit_length_set_bits(fmpz_mpoly_t A, slong len,
                              flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_is_fmpz_poly(const fmpz_mpoly_t A, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_get_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, flint_bitcnt_t Abits,
      const fmpz * Bcoeffs, slong Blen, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, const fmpz_poly_t B,
                                          slong v, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void tuple_print(fmpz * alpha, slong n);

FLINT_DLL void tuple_saturate(fmpz * alpha, slong n, slong m);

FLINT_DLL void tuple_next(fmpz * alpha, slong n);

/*****************************************************************************/

typedef struct {
    fmpz_t constant;
    fmpz_t constant_den;        /* should be one after normal operations */
    fmpz_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fmpz_mpoly_factor_struct;

typedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_init_set_ui(f->constant, 1);
    fmpz_init_set_ui(f->constant_den, 1);
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}

FLINT_DLL void fmpz_mpoly_factor_init2(fmpz_mpoly_factor_t f,
                                      slong alloc, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_realloc(fmpz_mpoly_factor_t f,
                                      slong alloc, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_fit_length(fmpz_mpoly_factor_t f,
                                        slong len, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f,
                                                   const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_FACTOR_INLINE
slong fmpz_mpoly_factor_length(const fmpz_mpoly_factor_t f,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return f->num;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c,
                      const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_set(c, f->constant);
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c,
                      const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_set(fmpq_numref(c), f->constant);
    fmpz_set(fmpq_denref(c), f->constant_den);
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_get_base(fmpz_mpoly_t p, const fmpz_mpoly_factor_t f,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mpoly_set(p, f->poly + i, ctx);
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_swap_base(fmpz_mpoly_t p, fmpz_mpoly_factor_t f,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mpoly_swap(p, f->poly + i, ctx);
}

FMPZ_MPOLY_FACTOR_INLINE
slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}

FLINT_DLL void fmpz_mpoly_factor_set(fmpz_mpoly_factor_t f,
                      const fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_cmp(const fmpz_mpoly_factor_t f,
                      const fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_print_pretty(const fmpz_mpoly_factor_t f,
                               const char ** vars, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_content(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_factor_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_factor_struct t = *f;
   *f = *g;
   *g = t;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_set_fmpz(fmpz_mpoly_factor_t f, const fmpz_t a,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    f->num = 0;
    fmpz_set(f->constant, a);
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_zero(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    f->num = 0;
    fmpz_zero(f->constant);
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_one(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    f->num = 0;
    fmpz_one(f->constant);
}

FLINT_DLL void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_expand(fmpz_mpoly_t A,
                      const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_bound_si(fmpz_t B, const fmpz_t A,
                                              const slong * degs, slong nvars);

FMPZ_MPOLY_FACTOR_INLINE
int fmpz_mpoly_factor_matches(const fmpz_mpoly_t A,
                       const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    int matches;
    fmpz_mpoly_t T;
    fmpz_mpoly_init(T, ctx);
    matches = fmpz_mpoly_factor_expand(T, f, ctx);
    matches = matches && fmpz_mpoly_equal(T, A, ctx);
    fmpz_mpoly_clear(T, ctx);
    return matches;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_append_fmpz_swap(fmpz_mpoly_factor_t f,
                    fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mpoly_swap(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_factor_append_ui(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A,
                                           ulong e, const fmpz_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mpoly_set(f->poly + i, A, ctx);
    fmpz_set_ui(f->exp + i, e);
    f->num = i + 1;
}

/*****************************************************************************/

typedef struct
{
    fmpz_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mpolyv_struct;

typedef fmpz_mpolyv_struct fmpz_mpolyv_t[1];

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpolyv_init(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpolyv_swap(fmpz_mpolyv_t A, fmpz_mpolyv_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void fmpz_mpolyv_clear(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyv_print_pretty(const fmpz_mpolyv_t poly,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyv_fit_length(fmpz_mpolyv_t A, slong length,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyv_set_coeff(fmpz_mpolyv_t A, slong i,
                                   fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_to_mpolyv(fmpz_mpolyv_t A, const fmpz_mpoly_t B,
                        const fmpz_mpoly_t xalpha, const fmpz_mpoly_ctx_t ctx);


FLINT_DLL void fmpz_mpoly_from_mpolyv(fmpz_mpoly_t A, flint_bitcnt_t Abits,
                            const fmpz_mpolyv_t B, const fmpz_mpoly_t xalpha,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_vec_content_mpoly(fmpz_mpoly_t g,
          const fmpz_mpoly_struct * A, slong Alen, const fmpz_mpoly_ctx_t ctx);


/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    fmpq_poly_struct * inv_prod_dbetas;
    fmpq_poly_struct * dbetas;
    fmpz_mpoly_struct * prod_mbetas;
    fmpz_mpolyv_struct * prod_mbetas_coeffs;
    fmpz_mpoly_struct * mbetas;
    fmpz_mpoly_struct * deltas;
    fmpz_mpoly_struct * xalpha;
    fmpz_mpoly_struct * q;
    fmpz_mpoly_univar_struct * U;
    fmpz_mpoly_geobucket_struct * G;
    fmpz_mpoly_struct * qt;
    fmpz_mpoly_struct * newt;
    fmpz_mpolyv_struct * delta_coeffs;
    fmpq_poly_t dtq, S, R;
} fmpz_mpoly_pfrac_struct;

typedef fmpz_mpoly_pfrac_struct fmpz_mpoly_pfrac_t[1];

FLINT_DLL int fmpz_mpoly_pfrac_init(fmpz_mpoly_pfrac_t I, flint_bitcnt_t bits,
                            slong r, slong w, const fmpz_mpoly_struct * betas,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_pfrac_clear(fmpz_mpoly_pfrac_t I,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_pfrac(slong l, fmpz_mpoly_t t, const slong * degs,
                             fmpz_mpoly_pfrac_t I, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_hlift(slong m, fmpz_mpoly_struct * f, slong r,
                const fmpz * alpha, const fmpz_mpoly_t A, const slong * degs,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_get_lead0(fmpz_mpoly_t c, const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set_lead0(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fmpz_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_bpoly_struct;

typedef fmpz_bpoly_struct fmpz_bpoly_t[1];

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_bpoly_init(fmpz_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_bpoly_swap(fmpz_bpoly_t A, fmpz_bpoly_t B)
{
    fmpz_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fmpz_bpoly_clear(fmpz_bpoly_t A);

FLINT_DLL void fmpz_bpoly_realloc(fmpz_bpoly_t A, slong len);

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_bpoly_fit_length(fmpz_bpoly_t A, slong len)
{
    if (A->alloc < len)
        fmpz_bpoly_realloc(A, len);
}

FLINT_DLL void fmpz_bpoly_print_pretty(fmpz_bpoly_t A,
                                         const char * var0, const char * var1);

FMPZ_MPOLY_FACTOR_INLINE
fmpz_poly_struct * fmpz_bpoly_lead(fmpz_bpoly_t A)
{
    return A->coeffs + A->length - 1;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_bpoly_zero(fmpz_bpoly_t A)
{
    A->length = 0;
}

FMPZ_MPOLY_FACTOR_INLINE
slong fmpz_bpoly_degree0(const fmpz_bpoly_t A)
{
    return A->length - 1;
}

FLINT_DLL slong fmpz_bpoly_degree1(const fmpz_bpoly_t A);

FLINT_DLL void fmpz_bpoly_set_coeff(fmpz_bpoly_t A, slong exp0, slong exp1,
                                                               const fmpz_t c);

FLINT_DLL void fmpz_mpoly_set_fmpz_bpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_bpoly_t B,
    slong var0,
    slong var1,
    const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_bpoly(
    fmpz_bpoly_t A,
    const fmpz_mpoly_t B,
    slong var0,
    slong var1,
    const fmpz_mpoly_ctx_t ctx);

typedef struct
{
    fmpz_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_tpoly_struct;

typedef fmpz_tpoly_struct fmpz_tpoly_t[1];

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_tpoly_init(fmpz_tpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_tpoly_swap(fmpz_tpoly_t A, fmpz_tpoly_t B)
{
    fmpz_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fmpz_tpoly_fit_length(fmpz_tpoly_t A, slong len);

FLINT_DLL void fmpz_tpoly_clear(fmpz_tpoly_t A);

FLINT_DLL void fmpz_bpoly_factor(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B);

FLINT_DLL int fmpz_bpoly_factor_ordered(
    fmpz_poly_t c,
    fmpz_tpoly_t F,
    fmpz_bpoly_t B,
    const fmpz_t alpha,
    const fmpz_poly_factor_t Bevalf);

/*****************************************************************************/

FMPZ_MPOLY_FACTOR_INLINE
void fmpz_mpoly_unit_normalize(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_mpoly_neg(A, A, ctx);
}

FLINT_DLL int _fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f,
                   fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx);


FLINT_DLL int fmpz_mpoly_factor_lcc_wang(fmpz_mpoly_struct * lc_divs,
                        const fmpz_mpoly_factor_t lcAfac, const fmpz_t Auc,
                    const fmpz_poly_struct * Auf, slong r, const fmpz * alpha,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_irred_zassenhaus(fmpz_mpolyv_t fac,
       const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx, zassenhaus_prune_t Z);

FLINT_DLL int fmpz_mpoly_factor_irred_wang(fmpz_mpolyv_t fac,
                      const fmpz_mpoly_t A, const fmpz_mpoly_factor_t lcAfac,
          int lcAfac_irred, const fmpz_mpoly_t lcA, const fmpz_mpoly_ctx_t ctx,
                    flint_rand_t state, zassenhaus_prune_t Z, int allow_shift);

FLINT_DLL int fmpz_mpoly_factor_irred_zippel(fmpz_mpolyv_t fac,
                    const fmpz_mpoly_t A, const fmpz_mpoly_factor_t lcAfac,
          int lcAfac_irred, const fmpz_mpoly_t lcA, const fmpz_mpoly_ctx_t ctx,
                                     flint_rand_t state, zassenhaus_prune_t Z);

FLINT_DLL int fmpz_mpoly_factor_irred(fmpz_mpoly_factor_t f,
                                const fmpz_mpoly_ctx_t ctx, unsigned int algo);

FLINT_DLL int fmpz_mpoly_factor_zassenhaus(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_wang(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_zippel(fmpz_mpoly_factor_t f,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_evaluate_rest_fmpz(fmpz * E,
         slong * starts, slong * ends, slong * stops, ulong * es,
          const fmpz * Acoeffs, const ulong * Aexps, slong Alen, slong var,
            const fmpz * alphas, const slong * offsets, const slong * shifts,
                                             slong N, ulong mask, slong nvars);

FLINT_DLL void _fmpz_mpoly_eval_rest_to_poly(fmpz_poly_t E,
        const fmpz_mpoly_t A, const fmpz * alphas, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_lcc_kaltofen_step(
    fmpz_mpoly_struct * divs,   /* length r */
    slong r,
    fmpz_mpoly_factor_t Af, /* squarefree factorization of A */
    const fmpz_poly_struct * Au,
    slong v,                      /* minor bivar var*/
    const fmpz * alphas,
    const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_lcc_kaltofen(
    fmpz_mpoly_struct * divs,
    const fmpz_mpoly_factor_t lcAf_,
    const fmpz_mpoly_t A,
    slong r,
    const fmpz * alpha,
    slong * degs,
    const fmpz_poly_factor_t uf,
    const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_evaluate_rest_except_one(
    fmpz_poly_t e,
    const fmpz_mpoly_t A,
    const fmpz * alphas,
    slong v,
    const fmpz_mpoly_ctx_t ctx);

/****************************************************************************/

FLINT_DLL void fmpz_mpoly_compression_do(fmpz_mpoly_t L,
                     const fmpz_mpoly_ctx_t Lctx, fmpz * Acoeffs, slong Alen,
                                                        mpoly_compression_t M);

FLINT_DLL void fmpz_mpoly_compression_undo(fmpz_mpoly_t A, flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx, fmpz_mpoly_t L, const fmpz_mpoly_ctx_t Lctx,
                                                        mpoly_compression_t M);

#ifdef __cplusplus
}
#endif

#endif

