/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MPOLY_FACTOR_H
#define NMOD_MPOLY_FACTOR_H

#ifdef NMOD_MPOLY_FACTOR_INLINES_C
#define NMOD_MPOLY_FACTOR_INLINE
#else
#define NMOD_MPOLY_FACTOR_INLINE static inline
#endif

#include "nmod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif


void nmod_mpoly_get_bpoly(n_bpoly_t A, const nmod_mpoly_t B,
                           slong var0, slong var1, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_bpoly(nmod_mpoly_t A, flint_bitcnt_t Abits,
        const n_bpoly_t B, slong var0, slong var1, const nmod_mpoly_ctx_t ctx);

int n_bpoly_mod_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    nmod_t ctx);

void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx);

/*****************************************************************************/

int nmod_mat_is_reduced(const nmod_mat_t N);

void nmod_mat_init_nullspace_tr(nmod_mat_t X, nmod_mat_t tmp);

/*****************************************************************************/

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
	f->constant = 1;
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}

void nmod_mpoly_factor_init2(nmod_mpoly_factor_t f, slong alloc,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_realloc(nmod_mpoly_factor_t f, slong alloc,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_fit_length(nmod_mpoly_factor_t f, slong len,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
slong nmod_mpoly_factor_length(const nmod_mpoly_factor_t f,
                                                    const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    return f->num;
}

NMOD_MPOLY_FACTOR_INLINE
ulong nmod_mpoly_factor_get_constant_ui(const nmod_mpoly_factor_t f,
                                                    const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    return f->constant;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_get_base(nmod_mpoly_t p, const nmod_mpoly_factor_t f,
                                           slong i, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    nmod_mpoly_set(p, f->poly + i, ctx);
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_swap_base(nmod_mpoly_t p, nmod_mpoly_factor_t f,
                                           slong i, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    nmod_mpoly_swap(p, f->poly + i, ctx);
}

slong nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t FLINT_UNUSED(ctx));

void nmod_mpoly_factor_append_ui(nmod_mpoly_factor_t f,
                    const nmod_mpoly_t A, ulong e, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_append_fmpz(nmod_mpoly_factor_t f,
             const nmod_mpoly_t A, const fmpz_t e, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_set(nmod_mpoly_factor_t f,
                      const nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_print_pretty(const nmod_mpoly_factor_t f,
                               const char ** vars, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_content(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_separable(nmod_mpoly_factor_t f,
                    const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, int sep);

int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_cmp(
    const nmod_mpoly_factor_t A,
    const nmod_mpoly_factor_t B,
    const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_expand(nmod_mpoly_t A,
                      const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
int nmod_mpoly_factor_matches(const nmod_mpoly_t a, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
    int matches;
    nmod_mpoly_t t;
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_factor_expand(t, f, ctx);
    matches = nmod_mpoly_equal(t, a, ctx);
    nmod_mpoly_clear(t, ctx);
    return matches;
}

int nmod_mpoly_factor_fix_units(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g,
                                                    const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
   nmod_mpoly_factor_struct t = *f;
   *f = *g;
   *g = t;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_one(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
	f->constant = 1;
	f->num = 0;
}

void _nmod_mpoly_get_lead0(
    nmod_mpoly_t c,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_set_lead0(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t c,
    const nmod_mpoly_ctx_t ctx);

/* n_poly_vec ****************************************************************/

slong _n_poly_vec_max_degree(const n_poly_struct * A, slong Alen);

void _n_poly_vec_mul_nmod_intertible(n_poly_struct * A,
                                          slong Alen, ulong c, nmod_t ctx);

void _n_poly_vec_mod_mul_poly(n_poly_struct * A, slong Alen,
                                           const n_poly_t g, const nmod_t ctx);

void _n_poly_vec_mod_divexact_poly(n_poly_struct * A, slong Alen,
                                                 const n_poly_t g, nmod_t ctx);

void _n_poly_vec_mod_content(n_poly_t g, const n_poly_struct * A,
                                                       slong Alen, nmod_t ctx);

void _n_poly_vec_mod_remove_content(n_poly_t g, n_poly_struct * A,
                                                       slong Alen, nmod_t ctx);

/* polyun ********************************************************************/

void nmod_mpoly_get_polyu1n(n_polyun_t A, const nmod_mpoly_t B,
                          slong varx, slong vary, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_polyu1n(nmod_mpoly_t B, const n_polyun_t A,
                          slong varx, slong vary, const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} nmod_mpolyv_struct;

typedef nmod_mpolyv_struct nmod_mpolyv_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpolyv_init(nmod_mpolyv_t A, const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpolyv_swap(nmod_mpolyv_t A, nmod_mpolyv_t B,
                                                    const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
   nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyv_print_pretty(const nmod_mpolyv_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyv_fit_length(nmod_mpolyv_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyv_set_coeff(nmod_mpolyv_t A, slong i,
                                   nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_to_mpolyv(nmod_mpolyv_t A, const nmod_mpoly_t B,
                        const nmod_mpoly_t xalpha, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_from_mpolyv(nmod_mpoly_t A, flint_bitcnt_t Abits,
                            const nmod_mpolyv_t B, const nmod_mpoly_t xalpha,
                                                   const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_vec_content_mpoly(nmod_mpoly_t g,
          const nmod_mpoly_struct * A, slong Alen, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_vec_divexact_mpoly(nmod_mpoly_struct * A,
                 slong Alen, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_vec_mul_mpoly(nmod_mpoly_struct * A,
                 slong Alen, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

int _nmod_mpoly_factor_separable(nmod_mpoly_factor_t f,
                    const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, int sep);

int nmod_mpoly_factor_lcc_wang(nmod_mpoly_struct * lc_divs,
             const nmod_mpoly_factor_t lcAfac, const n_poly_t Auc,
             const n_bpoly_struct * Auf, slong r, const n_poly_struct * alpha,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_irred_smprime_zassenhaus(nmod_mpolyv_t fac,
         const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_medprime_zassenhaus(nmod_mpolyv_t fac,
         const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_lgprime_zassenhaus(nmod_mpolyv_t fac,
         const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_smprime_wang(nmod_mpolyv_t fac,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_medprime_wang(nmod_mpolyv_t Af,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_lgprime_wang(nmod_mpolyv_t Af,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_smprime_zippel(nmod_mpolyv_t fac,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_medprime_zippel(nmod_mpolyv_t Af,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_irred_lgprime_zippel(nmod_mpolyv_t Af,
       const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac,
       const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state);

/*****************************************************************************/

void nmod_mpoly_compression_do(nmod_mpoly_t L,
                 const nmod_mpoly_ctx_t Lctx, ulong * Acoeffs, slong Alen,
                                                        mpoly_compression_t M);

void nmod_mpoly_compression_undo(nmod_mpoly_t A,
             flint_bitcnt_t Abits, const nmod_mpoly_ctx_t Actx, nmod_mpoly_t L,
                           const nmod_mpoly_ctx_t Lctx, mpoly_compression_t M);

/*****************************************************************************/

int nmod_mpolyu_is_canonical(const nmod_mpolyu_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu3_print_pretty(const nmod_mpolyu_t A,
                    const char * var0, const char * var1, const char * var2,
                               const char ** vars, const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    n_poly_struct * inv_prod_dbetas;
    nmod_mpoly_struct * inv_prod_dbetas_mvar;
    n_poly_struct * dbetas;
    nmod_mpoly_struct * dbetas_mvar;
    nmod_mpoly_struct * prod_mbetas;
    nmod_mpolyv_struct * prod_mbetas_coeffs;
    nmod_mpoly_struct * mbetas;
    nmod_mpoly_struct * deltas;
    nmod_mpoly_struct * xalpha;
    nmod_mpoly_struct * q;
    nmod_mpoly_geobucket_struct * G;
    nmod_mpoly_struct * qt;
    nmod_mpoly_struct * newt;
    nmod_mpolyv_struct * delta_coeffs;
    nmod_mpoly_t T;
    nmod_mpoly_t Q;
    nmod_mpoly_t R;
} nmod_mpoly_pfrac_struct;

typedef nmod_mpoly_pfrac_struct nmod_mpoly_pfrac_t[1];


int nmod_mpoly_pfrac_init(nmod_mpoly_pfrac_t Iv, flint_bitcnt_t bits,
                         slong l, slong r, const nmod_mpoly_struct * betas,
                          const ulong * alpha, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_pfrac_clear(nmod_mpoly_pfrac_t Iv,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_pfrac(slong r, nmod_mpoly_t t, const slong * deg,
                             nmod_mpoly_pfrac_t Iv, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_hlift(slong m, nmod_mpoly_struct * f, slong r,
            const ulong * alpha, const nmod_mpoly_t A, const slong * degs,
                                                   const nmod_mpoly_ctx_t ctx);

int n_bpoly_mod_pfrac(slong r, n_bpoly_struct * C,
            slong * C_deg1_bound, n_bpoly_t A, n_bpoly_struct * B, nmod_t mod);

int n_bpoly_mod_hlift2(n_bpoly_t A, n_bpoly_t B0, n_bpoly_t B1,
                              ulong alpha, slong degree_inner, nmod_t mod,
                                                      n_poly_bpoly_stack_t St);

int n_bpoly_mod_hlift2_cubic(n_bpoly_t A, n_bpoly_t B0, n_bpoly_t B1,
                               ulong alpha, slong degree_inner, nmod_t ctx,
                                nmod_eval_interp_t E, n_poly_bpoly_stack_t St);

int n_bpoly_mod_hlift(slong r, n_bpoly_t A, n_bpoly_struct * B,
                              ulong alpha, slong degree_inner, nmod_t mod,
                                                      n_poly_bpoly_stack_t St);

int n_bpoly_mod_hlift_cubic(slong r, n_bpoly_t A, n_bpoly_struct * B,
                               ulong alpha, slong degree_inner, nmod_t mod,
                                nmod_eval_interp_t E, n_poly_bpoly_stack_t St);

int n_polyu3_mod_hlift(slong r, n_polyun_struct * BB,  n_polyu_t A,
           n_polyu_struct * B, ulong beta, slong degree_inner, nmod_t ctx);

int nmod_mpoly_hlift_zippel(slong m, nmod_mpoly_struct * B, slong r,
            const ulong * alpha, const nmod_mpoly_t A, const slong * degs,
                               const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpoly_factor_algo(nmod_mpoly_factor_t f,
          const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, unsigned int algo);

int nmod_mpoly_factor_zassenhaus(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_wang(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_factor_zippel(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_evaluate_rest_n_poly(n_poly_struct * E,
    slong * starts, slong * ends, slong * stops, ulong * es,
    const ulong * Acoeffs, const ulong * Aexps, slong Alen, slong var,
    const n_poly_struct * alphas, const slong * offsets, const slong * shifts,
                                 slong N, ulong mask, slong nvars, nmod_t ctx);

void _nmod_mpoly_eval_rest_to_n_bpoly(n_bpoly_t E,
                      const nmod_mpoly_t A, const n_poly_struct * alphabetas,
                                                   const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_set_n_bpoly_var1_zero(nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const n_bpoly_t B, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

/* gcd ***********************************************************************/

int nmod_mpolyl_gcdp_zippel_smprime(nmod_mpoly_t G, nmod_mpoly_t Abar,
                nmod_mpoly_t Bbar, nmod_mpoly_t A, nmod_mpoly_t B, slong var,
                               const nmod_mpoly_ctx_t ctx, flint_rand_t state);

int nmod_mpolyl_gcds_zippel(nmod_mpoly_t G, const ulong * Gmarks,
        slong Gmarkslen, nmod_mpoly_t A, nmod_mpoly_t B, slong *perm,
        slong l, slong var, const nmod_mpoly_ctx_t ctx,  flint_rand_t state,
                          slong * Gdegbound, n_poly_t Amarks, n_poly_t Bmarks);

/* zip helpers ***************************************************************/

void mpoly_monomial_evals_nmod(n_poly_t EH, const ulong * Aexps,
         slong Alen, flint_bitcnt_t Abits, n_poly_struct * alpha_caches,
         slong start, slong stop, const mpoly_ctx_t mctx, const nmod_t fpctx);

void mpoly1_monomial_evals_nmod(n_polyun_t EH, const ulong * Aexps,
                flint_bitcnt_t Abits, const ulong * Amarks, slong Amarkslen,
                n_poly_struct * alpha_caches, slong m,  const mpoly_ctx_t mctx,
                                                           const nmod_t fpctx);

void mpoly2_monomial_evals_nmod(n_polyun_t EH, const ulong * Aexps,
                flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen,
                n_poly_struct * alpha_caches, slong m, const mpoly_ctx_t mctx,
                                                           const nmod_t fpctx);

void n_polyun_zip_start(n_polyun_t Z, n_polyun_t H, slong req_images);

int n_polyu2n_add_zip_must_match(n_polyun_t Z, const n_bpoly_t A,
                                                             slong cur_length);

int n_polyun_zip_solve(nmod_mpoly_t A, n_polyun_t Z, n_polyun_t H,
                                     n_polyun_t M, const nmod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
