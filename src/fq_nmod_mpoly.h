/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_H
#define FQ_NMOD_MPOLY_H

#ifdef FQ_NMOD_MPOLY_INLINES_C
#define FQ_NMOD_MPOLY_INLINE
#else
#define FQ_NMOD_MPOLY_INLINE static inline
#endif

#include "nmod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

/*
    context object for fq_nmod_mpoly
*/
typedef struct
{
    mpoly_ctx_t minfo;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpoly_ctx_struct;

typedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1];

/*
    fq_nmod_mpoly_t
    sparse multivariates with fq_nmod coefficients
*/
typedef struct {
    mp_limb_t * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} fq_nmod_mpoly_struct;

typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];


/* Internal type definitions *************************************************/

/*
    fq_nmod_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fq_nmod_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fq_nmod_mpoly_univar_struct;

typedef fq_nmod_mpoly_univar_struct fq_nmod_mpoly_univar_t[1];

/*
    fq_nmod_mpolyu_t
    sparse univariates with fq_nmod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   fq_nmod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} fq_nmod_mpolyu_struct;

typedef fq_nmod_mpolyu_struct fq_nmod_mpolyu_t[1];

/*
    fq_nmod_mpolyn_t
    multivariates with fq_nmod_poly_t coefficients
*/
typedef struct
{
   n_fq_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fq_nmod_mpolyn_struct;
typedef fq_nmod_mpolyn_struct fq_nmod_mpolyn_t[1];

/*
    fq_nmod_mpolyun_t
    sparse univariates with fq_nmod_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    fq_nmod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;   /* default bits to construct coeffs */
} fq_nmod_mpolyun_struct;

typedef fq_nmod_mpolyun_struct fq_nmod_mpolyun_t[1];

/* Embeddings ****************************************************************/

/* see fq_nmod_mpoly/fq_nmod_embed.c for more info */

typedef struct bad_fq_nmod_embed
{
    const fq_nmod_ctx_struct * smctx; /* modulus is f */
    fq_nmod_poly_t phi_sm;      /* phi as an element of F_p[theta][x] */
    fq_nmod_poly_t h;
    n_fq_poly_t h_as_n_fq_poly;
    const fq_nmod_ctx_struct * lgctx; /* modulus is g */
    fq_nmod_t theta_lg;         /* theta as an element of F_p[phi]/g(phi) */
    fq_nmod_t x_lg;             /* x as an element of F_p[phi]/g(phi) */
    nmod_mat_t lg_to_sm_mat;
    nmod_mat_t sm_to_lg_mat;
} bad_fq_nmod_embed_struct;

typedef bad_fq_nmod_embed_struct bad_fq_nmod_embed_t[1];


void bad_fq_nmod_embed_clear(bad_fq_nmod_embed_t emb);

void bad_fq_nmod_embed_array_init(bad_fq_nmod_embed_struct * emb,
                     const fq_nmod_ctx_t bigctx, const fq_nmod_ctx_t smallctx);

void bad_fq_nmod_embed_sm_to_lg(fq_nmod_t out,
                       const fq_nmod_poly_t in, const bad_fq_nmod_embed_t emb);

void bad_fq_nmod_embed_lg_to_sm(fq_nmod_poly_t out,
                            const fq_nmod_t in, const bad_fq_nmod_embed_t emb);

void bad_n_fq_embed_sm_to_lg(mp_limb_t * out_, const n_poly_t in_,
                                                const bad_fq_nmod_embed_t emb);

void bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(fq_nmod_t out,
                            const n_poly_t in_, const bad_fq_nmod_embed_t emb);

void bad_n_fq_embed_lg_to_sm(n_poly_t out_, const mp_limb_t * in_,
                                                const bad_fq_nmod_embed_t emb);

void bad_fq_nmod_embed_fq_nmod_lg_to_n_fq_sm(n_poly_t out_,
                            const fq_nmod_t in, const bad_fq_nmod_embed_t emb);

void bad_n_fq_embed_sm_elem_to_lg(mp_limb_t * out,
                          const mp_limb_t * in, const bad_fq_nmod_embed_t emb);

void bad_fq_nmod_embed_sm_elem_to_lg(fq_nmod_t out,
                            const fq_nmod_t in, const bad_fq_nmod_embed_t emb);

typedef struct bad_fq_nmod_mpoly_embed_chooser
{
    bad_fq_nmod_embed_struct * embed;
    slong m; /* degree of the extension F_q / F_p */
    slong n; /* degree of the extension F_q^n / F_q */
    slong k; /* index of current in embed */
    mp_limb_t p;
} bad_fq_nmod_mpoly_embed_chooser_struct;

typedef bad_fq_nmod_mpoly_embed_chooser_struct bad_fq_nmod_mpoly_embed_chooser_t[1];

bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_init(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);

void bad_fq_nmod_mpoly_embed_chooser_clear(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);

bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_next(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);


/* Context object ************************************************************/

void fq_nmod_mpoly_ctx_init_deg(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                 const ordering_t ord, mp_limb_t p, slong deg);

void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                              const ordering_t ord, const fq_nmod_ctx_t fqctx);

void fq_nmod_mpoly_ctx_init_rand(fq_nmod_mpoly_ctx_t ctx,
                                       flint_rand_t state, slong max_nvars,
                                          flint_bitcnt_t p_bits, slong deg_bound);

void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_ctx_nvars(const fq_nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FQ_NMOD_MPOLY_INLINE
ordering_t fq_nmod_mpoly_ctx_ord(const fq_nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}


/*  Memory management ********************************************************/

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->coeffs_alloc > 0)
        flint_free(A->coeffs);

    if (A->exps_alloc > 0)
        flint_free(A->exps);
}

void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc,
                           flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_fit_length_fit_bits(fq_nmod_mpoly_t A,
                slong len, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_fit_length_reset_bits(fq_nmod_mpoly_t A,
                slong len, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void _fq_nmod_mpoly_fit_length(
    mp_limb_t ** coeffs,
    slong * coeffs_alloc,
    slong d,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length)
{
    if (d*length > *coeffs_alloc)
    {
        *coeffs_alloc = FLINT_MAX(d*length, *coeffs_alloc*2);
        *coeffs = (mp_limb_t *) flint_realloc(*coeffs,
                                              *coeffs_alloc*sizeof(mp_limb_t));
    }

    if (N*length > *exps_alloc)
    {
        *exps_alloc = FLINT_MAX(N*length, *exps_alloc*2);
        *exps = (mp_limb_t *) flint_realloc(*exps, *exps_alloc*sizeof(ulong));
    }
}

FQ_NMOD_MPOLY_INLINE
void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(fq_nmod_ctx_degree(ctx->fqctx)*newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, ctx->minfo)*newlen <= A->exps_alloc);
    A->length = newlen;
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_truncate(fq_nmod_mpoly_t A, slong newlen,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        A->length = newlen;
    }
}


/* Input/output **************************************************************/

int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fq_nmod_mpoly_fprint_pretty(FILE * file, const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx);
#endif

int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx);

/*  Basic manipulation *******************************************************/

void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_equal(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fq_nmod_mpoly_struct, *A, *B);
}


/* Constants *****************************************************************/

FQ_NMOD_MPOLY_INLINE
mp_limb_t * fq_nmod_mpoly_get_nonzero_n_fq(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == 1);
    FLINT_ASSERT(mpoly_monomial_is_zero(A->exps,
                                    mpoly_words_per_exp(A->bits, ctx->minfo)));
    return A->coeffs;
}

int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_n_fq(fq_nmod_mpoly_t A,
                           const mp_limb_t * c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_fmpz(fq_nmod_mpoly_t A, const fmpz_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   _fq_nmod_mpoly_set_length(A, 0, ctx);
}

void fq_nmod_mpoly_one(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* Degrees *******************************************************************/

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degrees_si(slong * degs, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_total_degree_fmpz(fmpz_t td, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_used_vars(int * used, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}

/* Coefficients **************************************************************/

void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c,
                          const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A,
                                  const fq_nmod_t c, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c,
                             const fq_nmod_mpoly_t A, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c,
                                const fq_nmod_mpoly_t A, const ulong * exp,
                                               const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                       const fq_nmod_t c, const fmpz * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                     const fq_nmod_t c, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A,
                                      const fq_nmod_t c, const ulong * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C,
              const fq_nmod_mpoly_t A, const slong * vars, const ulong * exps,
                                  slong length, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE mp_limb_t * _fq_nmod_mpoly_leadcoeff(
                        const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

FQ_NMOD_MPOLY_INLINE int fq_nmod_mpoly_is_monic(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length > 0 &&
           _n_fq_is_one(A->coeffs + 0, fq_nmod_ctx_degree(ctx->fqctx));
}

/* conversion ****************************************************************/

int fq_nmod_mpoly_is_fq_nmod_poly(const fq_nmod_mpoly_t A,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_get_fq_nmod_poly(fq_nmod_poly_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_set_fq_nmod_poly(fq_nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const fq_nmod_struct * Bcoeffs,
                         slong Blen, slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_fq_nmod_poly(fq_nmod_mpoly_t A,
             const fq_nmod_poly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

/* comparison ****************************************************************/

int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, slong new_length,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_term_coeff_fq_nmod(fq_nmod_mpoly_t A,
                    slong i, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_term_exp_ui(ulong * exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_term_exp_si(slong * exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

ulong fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A,
                            slong i, slong var, const fq_nmod_mpoly_ctx_t ctx);

slong fq_nmod_mpoly_get_term_var_exp_si(const fq_nmod_mpoly_t A,
                            slong i, slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A,
                   slong i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A,
                    slong i, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A,
         const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_push_term_fq_nmod_ffmpz(fq_nmod_mpoly_t A,
         const fq_nmod_t c, const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A,
          const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_assert_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_radix_sort1(fq_nmod_mpoly_t A, slong left,
     slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask, slong d);

void _fq_nmod_mpoly_radix_sort(fq_nmod_mpoly_t A, slong left,
           slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask, slong d);

void _fq_nmod_mpoly_push_exp_ffmpz(fq_nmod_mpoly_t A,
                              const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_push_exp_pfmpz(fq_nmod_mpoly_t A,
                            fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_push_exp_ui(fq_nmod_mpoly_t A,
                             const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state,
                 slong length, ulong exp_bound, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state,
              slong length, ulong * exp_bounds, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state,
            slong length, flint_bitcnt_t exp_bits, const fq_nmod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

slong _fq_nmod_mpoly_add(
                         mp_limb_t * coeff1,       ulong * exp1,
                         mp_limb_t * coeff2, const ulong * exp2, slong len2,
                         mp_limb_t * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);

void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_add_n_fq(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const mp_limb_t * c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_add(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_scalar_mul_n_fq(fq_nmod_mpoly_t A,
  const fq_nmod_mpoly_t B, const mp_limb_t * c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_scalar_addmul_fq_nmod(fq_nmod_mpoly_t A,
        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_t e,
                                                const fq_nmod_mpoly_ctx_t ctx);

/* Differentiation **********************************************************/

void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

void fq_nmod_mpoly_evaluate_one_fq_nmod(fq_nmod_mpoly_t A,
                    const fq_nmod_mpoly_t B, slong var, const fq_nmod_t val,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_eval_all_fq_nmod(fq_nmod_t ev,
                  const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
                       flint_bitcnt_t Abits, fq_nmod_struct * const * alphas,
                            const mpoly_ctx_t mctx, const fq_nmod_ctx_t fqctx);

void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, const fq_nmod_mpoly_t A,
                 fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A,
                     const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_compose_mat(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fmpz_mat_t M,
              const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC);

int fq_nmod_mpoly_compose_fq_nmod_mpoly_geobucket(fq_nmod_mpoly_t A,
              const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C,
              const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC);

int fq_nmod_mpoly_compose_fq_nmod_mpoly_horner(fq_nmod_mpoly_t A,
              const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C,
              const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC);

int fq_nmod_mpoly_compose_fq_nmod_mpoly(fq_nmod_mpoly_t A,
              const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C,
              const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC);

void fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(fq_nmod_mpoly_t A,
                             const fq_nmod_mpoly_t B, const slong * c,
              const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC);


/* Multiplication ************************************************************/

void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                       const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t poly1,
                    const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_mul_johnson(
    fq_nmod_mpoly_t A,
    const mp_limb_t * Bcoeffs,
    const ulong * Bexps,
    slong Blen,
    const mp_limb_t * Ccoeffs,
    const ulong * Cexps,
    slong Clen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const fq_nmod_ctx_t ctx);


/* Powering ******************************************************************/

int fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);


/* Division ******************************************************************/

int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q,
                         const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q,
                                  fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A,
                                  fq_nmod_mpoly_struct * const * B, slong len,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_divexact(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_divides(Q, A, B, ctx))
        return;

    flint_throw(FLINT_ERROR, "fq_nmod_mpoly_divexact: nonexact division");
}


int fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t poly1,
                  const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_div_monagan_pearce(fq_nmod_mpoly_t q,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_divrem_monagan_pearce(fq_nmod_mpoly_t q, fq_nmod_mpoly_t r,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_divrem_ideal_monagan_pearce(
                        fq_nmod_mpoly_struct ** q, fq_nmod_mpoly_t r,
            const fq_nmod_mpoly_t poly2, fq_nmod_mpoly_struct * const * poly3,
                                      slong len, const fq_nmod_mpoly_ctx_t ctx);

int _fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t A,
             const mp_limb_t * coeff2, const ulong * exp2, slong len2,
             const mp_limb_t * coeff3, const ulong * exp3, slong len3,
  flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);


/* Square root ***************************************************************/

int fq_nmod_mpoly_sqrt_heap(fq_nmod_mpoly_t Q,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_sqrt(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_sqrt_heap(Q, A, ctx);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_is_square(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int res;
    fq_nmod_mpoly_t Q;
    fq_nmod_mpoly_init(Q, ctx);
    res = fq_nmod_mpoly_sqrt_heap(Q, A, ctx);
    fq_nmod_mpoly_clear(Q, ctx);
    return res;
}

int fq_nmod_mpoly_quadratic_root(fq_nmod_mpoly_t Q,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

void fq_nmod_mpoly_term_content(fq_nmod_mpoly_t M,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_content_vars(fq_nmod_mpoly_t g,
                     const fq_nmod_mpoly_t A, slong * vars, slong vars_length,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

int _fq_nmod_mpoly_gcd_algo(fq_nmod_mpoly_t G, fq_nmod_mpoly_t Abar,
        fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                             const fq_nmod_mpoly_ctx_t ctx, unsigned int algo);

int fq_nmod_mpoly_gcd_cofactors(fq_nmod_mpoly_t G,
          fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_gcd_hensel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_gcd_zippel2(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_deflation(fmpz * shift, fmpz * stride,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_deflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_inflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
       const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx);


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

void mpoly_void_ring_init_fq_nmod_mpoly_ctx(mpoly_void_ring_t R,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyl_lead_coeff(fq_nmod_mpoly_t c,
                                     const fq_nmod_mpoly_t A, slong num_vars,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_mpolyl_perm_deflate(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t lctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

void fq_nmod_mpoly_from_mpolyl_perm_inflate(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

int fq_nmod_mpolyl_gcd_zippel_smprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyl_gcd_zippel_lgprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyl_gcd_hensel_smprime(
    fq_nmod_mpoly_t G, slong Gdeg,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyl_content(fq_nmod_mpoly_t g,
                                    const fq_nmod_mpoly_t A, slong num_vars,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_pow_rmul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_to_fq_nmod_poly_deflate(fq_nmod_poly_t A,
                   const fq_nmod_mpoly_t B, slong var, const ulong * Bshift,
                         const ulong * Bstride, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_from_fq_nmod_poly_inflate(fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits, const fq_nmod_poly_t B, slong var, const ulong * Ashift,
                         const ulong * Astride, const fq_nmod_mpoly_ctx_t ctx);


int fq_nmod_mpoly_repack_bits(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                          flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_repack_bits_inplace(fq_nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx);


void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx,
                                                                    slong deg);


/* Univariates ***************************************************************/

void fq_nmod_mpoly_univar_init(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_univar_fit_length(fq_nmod_mpoly_univar_t A,
                                  slong length, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_univar_print_pretty(const fq_nmod_mpoly_univar_t A,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_univar_assert_canonical(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_zero(fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

void fq_nmod_mpoly_univar_set_coeff_ui(fq_nmod_mpoly_univar_t A,
              ulong e, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, flint_bitcnt_t Abits,
     const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A,
     const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_swap(fq_nmod_mpoly_univar_t A,
                       fq_nmod_mpoly_univar_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fq_nmod_mpoly_univar_struct, *A, *B);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_univar_degree_fits_si(const fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_univar_length(const fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_univar_get_term_exp_si(fq_nmod_mpoly_univar_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_get_term_coeff(fq_nmod_mpoly_t c,
        const fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong)A->length);
    fq_nmod_mpoly_set(c, A->coeffs + i, ctx);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_swap_term_coeff(fq_nmod_mpoly_t c,
              fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong)A->length);
    fq_nmod_mpoly_swap(c, A->coeffs + i, ctx);
}

int fq_nmod_mpoly_univar_pseudo_gcd(fq_nmod_mpoly_univar_t Gx,
        const fq_nmod_mpoly_univar_t Ax, const fq_nmod_mpoly_univar_t Bx,
                                               const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_univar_resultant(fq_nmod_mpoly_t R,
        const fq_nmod_mpoly_univar_t Ax, const fq_nmod_mpoly_univar_t Bx,
                                               const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_univar_discriminant(fq_nmod_mpoly_t D,
             const fq_nmod_mpoly_univar_t Fx, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_resultant(fq_nmod_mpoly_t R,
                        const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                    slong var, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpoly_discriminant(fq_nmod_mpoly_t R,
          const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx);

/* mpolyu ********************************************************************/

int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t poly,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, flint_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpolyu_swap(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                               const fq_nmod_mpoly_ctx_t uctx)
{
    fq_nmod_mpolyu_struct t = *B;
    *B = *A;
    *A = t;
}

void fq_nmod_mpolyu_zero(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

int fq_nmod_mpolyu_is_one(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length,
                                               const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_mpolyu_one(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_mpolyu_degrees_si(
    slong * degs,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_repack_bits_inplace(
    fq_nmod_mpolyu_t A,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s);

void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s);

int fq_nmod_mpolyu_content_mpoly(fq_nmod_mpoly_t g,
                      const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A,
                                   fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B,
                                               const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_mpolyu_evaluate_one_fq_nmod(fq_nmod_mpolyu_t E,
                        fq_nmod_mpolyu_t A, slong var, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A,
                                    ulong pow, const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_mpoly_to_mpolyu_perm_deflate(
                         fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx,
                 const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride);

void fq_nmod_mpoly_from_mpolyu_perm_inflate(
        fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx,
          const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx,
                const slong * perm, const ulong * shift, const ulong * stride);

int fq_nmod_mpolyuu_divides(fq_nmod_mpolyu_t Q,
          const fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B, slong nmainvars,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_divexact_mpoly_inplace(fq_nmod_mpolyu_t A,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_mul_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_mul_mpoly_inplace(fq_nmod_mpolyu_t A,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyu_gcdm_zippel(fq_nmod_mpolyu_t G,
             fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar,  fq_nmod_mpolyu_t A,
          fq_nmod_mpolyu_t B, fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate);

FQ_NMOD_MPOLY_INLINE mp_limb_t * fq_nmod_mpolyu_leadcoeff(
                       const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _fq_nmod_mpoly_leadcoeff(A->coeffs + 0, ctx);
}

/* mpolyn ********************************************************************/

void fq_nmod_mpolyn_init(fq_nmod_mpolyn_t A, flint_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_clear(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_swap(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B);

int fq_nmod_mpolyn_is_canonical(const fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_zero(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_one(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_is_zero(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_print_pretty(const fq_nmod_mpolyn_t A,
                             const char ** x_in, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_fit_length(fq_nmod_mpolyn_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_fit_bits(fq_nmod_mpolyn_t A, slong bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_set(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE mp_limb_t * fq_nmod_mpolyn_leadcoeff(
                             fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    n_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + d*(leadpoly->length - 1);
}

FQ_NMOD_MPOLY_INLINE n_poly_struct * fq_nmod_mpolyn_leadcoeff_poly(
                       const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* mpolyun *******************************************************************/

void fq_nmod_mpoly_to_mpolyn_perm_deflate(fq_nmod_mpolyn_t A,
                const fq_nmod_mpoly_ctx_t nctx,
                const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride);

void fq_nmod_mpoly_from_mpolyn_perm_inflate(fq_nmod_mpoly_t A,
                flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx,
                const fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t nctx,
                const slong * perm, const ulong * shift, const ulong * stride);

void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A,
                           flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyun_is_canonical(const fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B);

void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A,
                                  slong length, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_one(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_is_nonzero_fq_nmod(const fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyun_is_nonzero_fq_nmod(const fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_scalar_mul_fq_nmod(fq_nmod_mpolyn_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_scalar_mul_fq_nmod(fq_nmod_mpolyun_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s);

void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s);

void fq_nmod_mpolyn_mul_poly(fq_nmod_mpolyn_t A,
                        const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c,
                              const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t t);

void fq_nmod_mpolyun_mul_poly(fq_nmod_mpolyun_t A,
                         const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_divexact_poly(fq_nmod_mpolyn_t A,
                      const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c,
            const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t q, fq_nmod_poly_t r);

void fq_nmod_mpolyun_divexact_poly(fq_nmod_mpolyun_t A,
                        const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_content_poly(fq_nmod_poly_t g,
                            fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_content_poly(fq_nmod_poly_t g,
                           fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx);

slong fq_nmod_mpolyn_lastdeg(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_set(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyun_t B,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_cvtto_mpolyn(
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_cvtto_mpolyun(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyu_t B,
    slong k,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_cvtfrom_mpolyn(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyu_cvtfrom_mpolyun(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyun_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE n_poly_struct * fq_nmod_mpolyun_leadcoeff_poly(
                      const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fq_nmod_mpolyn_leadcoeff_poly(A->coeffs + 0, ctx);
}


/*** gcd *********************************************************************/

int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

void fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

nmod_gcds_ret_t fq_nmod_mpolyu_gcds_zippel(fq_nmod_mpolyu_t G,
                fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpolyu_t f,
                                    slong var, const fq_nmod_mpoly_ctx_t ctx,
                                     flint_rand_t randstate, slong * degbound);

int fq_nmod_mpolyu_gcdp_zippel_univar(fq_nmod_mpolyu_t G,
            fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A,
                            fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(fq_nmod_mpolyu_t G,
        fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyu_gcdp_zippel(fq_nmod_mpolyu_t G,
              fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A,
                  fq_nmod_mpolyu_t B, slong var, const fq_nmod_mpoly_ctx_t ctx,
                                                       flint_rand_t randstate);

int fq_nmod_mpolyn_gcd_brown_smprime(fq_nmod_mpolyn_t G,
         fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_gcd_brown_lgprime(fq_nmod_mpolyn_t G,
         fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_monomial_evals2_cache(n_fq_polyun_t E,
        const ulong * Aexps, flint_bitcnt_t Abits, slong Alen,
        const fq_nmod_struct * betas, slong m, const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_monomial_evals_cache(n_fq_poly_t E,
                    const ulong * Aexps, flint_bitcnt_t Abits, slong Alen,
                    const fq_nmod_struct * betas, slong start, slong stop,
                                                const fq_nmod_mpoly_ctx_t ctx);

void n_fq_bpoly_eval_step_sep(n_fq_bpoly_t E, n_fq_polyun_t cur,
                            const n_fq_polyun_t inc, const fq_nmod_mpoly_t A,
                                                      const fq_nmod_ctx_t ctx);

void n_fq_polyun_zip_start(n_fq_polyun_t Z, n_fq_polyun_t H,
                                    slong req_images, const fq_nmod_ctx_t ctx);

int n_fq_polyu2n_add_zip_must_match(n_fq_polyun_t Z,
              const n_fq_bpoly_t A, slong cur_length, const fq_nmod_ctx_t ctx);

int n_fq_polyun_zip_solve(fq_nmod_mpoly_t A, n_fq_polyun_t Z,
              n_fq_polyun_t H, n_fq_polyun_t M, const fq_nmod_mpoly_ctx_t ctx);

/* gcd_helper_eval_interp ****************************************************/

void nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E,
                            const fq_nmod_ctx_t fqctx, const nmod_mpolyn_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_,
                         nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx,
                            const fq_nmod_poly_t B, const fq_nmod_ctx_t fqctx);

int nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_,
     nmod_mpolyn_t F, nmod_mpolyn_t T, n_poly_t modulus,
     const nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A,  const fq_nmod_ctx_t fqctx);

void nmod_mpolyn_interp_lift_lg_bpoly(slong * lastdeg_,
              nmod_mpolyn_t F, const nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A,
                                              const fq_nmod_mpoly_ctx_t lgctx);

int nmod_mpolyn_interp_crt_lg_bpoly(slong * lastdeg, nmod_mpolyn_t F,
           nmod_mpolyn_t T, n_fq_poly_t modulus, const nmod_mpoly_ctx_t smctx,
                              n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t lgctx);

void nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E,
                      fq_nmod_mpoly_ctx_t ectx, nmod_mpolyn_t A, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg,
                nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx,
                           fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx);

int nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_,
            nmod_mpolyn_t F, nmod_mpolyn_t T, n_poly_t modulus, slong var,
            const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A,
                                               const fq_nmod_mpoly_ctx_t ectx);

void nmod_mpolyn_interp_reduce_lg_mpoly(fq_nmod_mpoly_t A,
                            nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ffctx,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A,
                          nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ffctx,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_lg_mpoly(nmod_mpolyn_t A,
                              const nmod_mpoly_ctx_t ctx, fq_nmod_mpoly_t Ap,
                                               const fq_nmod_mpoly_ctx_t ctxp);

void nmod_mpolyun_interp_lift_lg_mpolyu(nmod_mpolyun_t A,
                            const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t Ap,
                                               const fq_nmod_mpoly_ctx_t ctxp);

int nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg,
                          nmod_mpolyun_t F, nmod_mpolyun_t T, n_poly_t m,
                          const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                                              const fq_nmod_mpoly_ctx_t ffctx);

int nmod_mpolyn_interp_mcrt_lg_mpoly(slong * lastdeg_,
            nmod_mpolyn_t H, const nmod_mpoly_ctx_t smctx, const n_poly_t m,
                              const mp_limb_t * inv_m_eval,fq_nmod_mpoly_t A,
                                              const fq_nmod_mpoly_ctx_t lgctx);

int nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg,
               nmod_mpolyun_t H, const nmod_mpoly_ctx_t ctx, n_poly_t m,
                           fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctxp);

void fq_nmod_mpolyn_interp_reduce_sm_poly(fq_nmod_poly_t E,
                            const fq_nmod_mpolyn_t A, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_lift_sm_poly(fq_nmod_mpolyn_t A,
                        const fq_nmod_poly_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_interp_crt_sm_poly(slong * lastdeg_,
            fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, const fq_nmod_poly_t A,
            const fq_nmod_poly_t modulus, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_lift_sm_bpoly(fq_nmod_mpolyn_t F,
                                n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_interp_crt_sm_bpoly(slong * lastdeg,
                fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, const n_fq_bpoly_t A,
                const n_fq_poly_t modulus, n_fq_poly_t alphapow,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_reduce_sm_mpolyn(fq_nmod_mpolyn_t E,
                             fq_nmod_mpolyn_t A, slong var, fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_lift_sm_mpolyn(fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_interp_mcrt_sm_mpoly(slong * lastdeg_,
        fq_nmod_mpolyn_t F, fq_nmod_mpoly_t A, const n_fq_poly_t modulus,
                          n_fq_poly_t alphapow, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_interp_crt_sm_mpolyn(slong * lastdeg_,
        fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_mpolyn_t A,  slong var,
        fq_nmod_poly_t modulus, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E,
                        const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_,
      fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

int fq_nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_,
            fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus,
            const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyn_interp_lift_lg_bpoly(slong * lastdeg_,
         fq_nmod_mpolyn_t F, const fq_nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A,
               const fq_nmod_mpoly_ctx_t lgctx, const bad_fq_nmod_embed_t emb);

int fq_nmod_mpolyn_interp_crt_lg_bpoly(slong * lastdeg,
               fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, n_fq_poly_t modulus,
               const fq_nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A,
               const fq_nmod_mpoly_ctx_t lgctx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E,
               const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A, slong var,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg_,
                         fq_nmod_mpolyn_t A, slong var,
                         const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

int fq_nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_,
    fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus, slong var,
    const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A,
                  const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyun_interp_reduce_sm_mpolyu(fq_nmod_mpolyu_t B,
          fq_nmod_mpolyun_t A, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_lift_sm_mpoly(
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyun_interp_lift_sm_mpolyu(fq_nmod_mpolyun_t A,
                      const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyun_interp_crt_sm_mpolyu(slong * lastdeg,
            fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_mpolyu_t A,
       fq_nmod_poly_t modulus, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpolyn_interp_reduce_lg_mpoly(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A,
                        fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ectx,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyn_interp_lift_lg_mpoly(fq_nmod_mpolyn_t A,
                            const fq_nmod_mpoly_ctx_t ctx,fq_nmod_mpoly_t B,
                const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

void fq_nmod_mpolyun_interp_lift_lg_mpolyu(fq_nmod_mpolyun_t A,
                          const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

int fq_nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg,
                fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t m,
                const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

int fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg,
                        fq_nmod_mpolyun_t H, const fq_nmod_mpoly_ctx_t ctx,
                        fq_nmod_poly_t m, fq_nmod_mpolyu_t A,
                         const fq_nmod_mpoly_ctx_t ectx, bad_fq_nmod_embed_t emb);

/* geobuckets ****************************************************************/

typedef struct fq_nmod_mpoly_geobucket
{
    fq_nmod_mpoly_struct polys[FLINT_BITS/2];
    fq_nmod_mpoly_struct temps[FLINT_BITS/2];
    slong length;
} fq_nmod_mpoly_geobucket_struct;

typedef fq_nmod_mpoly_geobucket_struct fq_nmod_mpoly_geobucket_t[1];

void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p,
                   fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal consistency checks

******************************************************************************/

void fq_nmod_mpoly_remainder_strongtest(const fq_nmod_mpoly_t r, const fq_nmod_mpoly_t g, const fq_nmod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

