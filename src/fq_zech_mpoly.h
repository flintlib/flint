/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_MPOLY_H
#define FQ_ZECH_MPOLY_H

#ifdef FQ_ZECH_MPOLY_INLINES_C
#define FQ_ZECH_MPOLY_INLINE
#else
#define FQ_ZECH_MPOLY_INLINE static inline
#endif

#include "fq_zech_poly.h"
#include "fq_nmod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

FQ_ZECH_MPOLY_INLINE
nmod_t fq_zech_ctx_mod(const fq_zech_ctx_t ctx)
{
    return ctx->fq_nmod_ctx->mod;
}


/*  Type definitions *********************************************************/

/*
    context object for fq_zech_mpoly
*/
typedef struct
{
    mpoly_ctx_t minfo;
    fq_zech_ctx_t fqctx;
} fq_zech_mpoly_ctx_struct;

typedef fq_zech_mpoly_ctx_struct fq_zech_mpoly_ctx_t[1];

/*
    fq_zech_mpoly_t
    sparse multivariates with fq_zech coefficients
*/
typedef struct
{
    fq_zech_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;     /* number of bits per exponent */
} fq_zech_mpoly_struct;

typedef fq_zech_mpoly_struct fq_zech_mpoly_t[1];

/* Internal type definitions *************************************************/

/*
    fq_zech_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fq_zech_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fq_zech_mpoly_univar_struct;

typedef fq_zech_mpoly_univar_struct fq_zech_mpoly_univar_t[1];

/*
    fq_zech_mpolyu_t
    sparse univariates with fq_zech_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   fq_zech_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} fq_zech_mpolyu_struct;

typedef fq_zech_mpolyu_struct fq_zech_mpolyu_t[1];

/*
    fq_zech_mpolyn_t
    multivariates with fq_zech_poly_t coefficients
*/
typedef struct
{
   fq_zech_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fq_zech_mpolyn_struct;
typedef fq_zech_mpolyn_struct fq_zech_mpolyn_t[1];

/*
    fq_zech_mpolyun_t
    sparse univariates with fq_zech_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    fq_zech_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;   /* default bits to construct coeffs */
} fq_zech_mpolyun_struct;

typedef fq_zech_mpolyun_struct fq_zech_mpolyun_t[1];

/*
    fq_zech_mpoly_geobucket_t
    power of 4 increment
*/
typedef struct fq_zech_mpoly_geobucket
{
    fq_zech_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fq_zech_mpoly_geobucket_struct;

typedef fq_zech_mpoly_geobucket_struct fq_zech_mpoly_geobucket_t[1];


/* Context object ************************************************************/

void fq_zech_mpoly_ctx_init_deg(fq_zech_mpoly_ctx_t ctx, slong nvars,
                                 const ordering_t ord, mp_limb_t p, slong deg);

void fq_zech_mpoly_ctx_clear(fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_ctx_nvars(const fq_zech_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FQ_ZECH_MPOLY_INLINE
ordering_t fq_zech_mpoly_ctx_ord(const fq_zech_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}


/*  Memory management ********************************************************/

void fq_zech_mpoly_init(fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_init2(fq_zech_mpoly_t A, slong alloc,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_init3(fq_zech_mpoly_t A, slong alloc,
                              flint_bitcnt_t bits, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_realloc(fq_zech_mpoly_t A,
                                   slong alloc, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_fit_length(fq_zech_mpoly_t A, slong length,
                                                const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_fit_length(fq_zech_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N,
                                                    const fq_zech_ctx_t fqctx);

void fq_zech_mpoly_fit_length_reset_bits(fq_zech_mpoly_t A,
                slong len, flint_bitcnt_t bits, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_clear(fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
void _fq_zech_mpoly_set_length(fq_zech_mpoly_t A, slong newlen,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->alloc);
    A->length = newlen;
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_truncate(fq_zech_mpoly_t A, slong newlen,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        A->length = newlen;
    }
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_fit_bits(fq_zech_mpoly_t A, slong bits,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    if (A->bits < bits)
    {
        if (A->alloc != 0)
        {
            slong N = mpoly_words_per_exp(bits, ctx->minfo);
            ulong * t = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length,
                                                                   ctx->minfo);
            flint_free(A->exps);
            A->exps = t;
        }

        A->bits = bits;
    }
}


/* Input/output **************************************************************/

int fq_zech_mpoly_set_str_pretty(fq_zech_mpoly_t A, const char * str,
                               const char ** x, const fq_zech_mpoly_ctx_t ctx);

char * fq_zech_mpoly_get_str_pretty(const fq_zech_mpoly_t A,
                               const char ** x, const fq_zech_mpoly_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fq_zech_mpoly_fprint_pretty(FILE * file, const fq_zech_mpoly_t A, const char ** x, const fq_zech_mpoly_ctx_t ctx);
#endif

int fq_zech_mpoly_print_pretty(const fq_zech_mpoly_t A, const char ** x, const fq_zech_mpoly_ctx_t ctx);

/*  Basic manipulation *******************************************************/

void fq_zech_mpoly_gen(fq_zech_mpoly_t A, slong var,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_is_gen(const fq_zech_mpoly_t A,
                                     slong var, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_equal(const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_swap(fq_zech_mpoly_t A, fq_zech_mpoly_t B,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fq_zech_mpoly_struct, *A, *B);
}


/* Constants *****************************************************************/

int fq_zech_mpoly_is_fq_zech(const fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_fq_zech(fq_zech_t c, const fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_fq_zech(fq_zech_mpoly_t A,
                             const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_ui(fq_zech_mpoly_t A, ulong c,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_fq_zech_gen(fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_equal_fq_zech(const fq_zech_mpoly_t A,
                             const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_zero(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
   _fq_zech_mpoly_set_length(A, 0, ctx);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_one(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_mpoly_set_ui(A, 1, ctx);
}

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_is_zero(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

int fq_zech_mpoly_is_one(const fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);


/* Degrees *******************************************************************/

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_degrees_fit_si(const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_degrees_fmpz(fmpz ** degs, const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_degrees_si(slong * degs, const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_degree_fmpz(fmpz_t deg, const fq_zech_mpoly_t A, slong var,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_degree_si(const fq_zech_mpoly_t A, slong var,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_total_degree_fits_si(const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_total_degree_fmpz(fmpz_t td, const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_total_degree_si(const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}


/* Coefficients **************************************************************/

void fq_zech_mpoly_get_coeff_fq_zech_monomial(fq_zech_t c,
                          const fq_zech_mpoly_t A, const fq_zech_mpoly_t M,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_coeff_fq_zech_monomial(fq_zech_mpoly_t A,
                                  const fq_zech_t c, const fq_zech_mpoly_t M,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_coeff_fq_zech_fmpz(fq_zech_t c,
                             const fq_zech_mpoly_t A, fmpz * const * exp,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_coeff_fq_zech_ui(fq_zech_t c,
                                const fq_zech_mpoly_t A, const ulong * exp,
                                               const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_set_coeff_fq_zech_fmpz(fq_zech_mpoly_t A,
                                       const fq_zech_t c, const fmpz * exp,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_coeff_fq_zech_fmpz(fq_zech_mpoly_t A,
                                     const fq_zech_t c, fmpz * const * exp,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_coeff_fq_zech_ui(fq_zech_mpoly_t A,
                                      const fq_zech_t c, const ulong * exp,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_coeff_vars_ui(fq_zech_mpoly_t C,
              const fq_zech_mpoly_t A, const slong * vars, const ulong * exps,
                                  slong length, const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE fq_zech_struct * fq_zech_mpoly_leadcoeff(
                        const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}


/* comparison ****************************************************************/

int fq_zech_mpoly_cmp(const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);


/* container operations ******************************************************/

int fq_zech_mpoly_is_canonical(const fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_length(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    return A->length;
}

void fq_zech_mpoly_resize(fq_zech_mpoly_t A, slong new_length,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term_coeff_fq_zech(fq_zech_t c,
              const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_term_coeff_fq_zech(fq_zech_mpoly_t A,
                    slong i, const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_term_exp_fits_ui(const fq_zech_mpoly_t A, slong i,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_term_exp_fits_si(const fq_zech_mpoly_t A, slong i,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

void fq_zech_mpoly_get_term_exp_fmpz(fmpz ** exp,
              const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term_exp_ui(ulong * exp,
              const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term_exp_si(slong * exp,
              const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx);

ulong fq_zech_mpoly_get_term_var_exp_ui(const fq_zech_mpoly_t A,
                            slong i, slong var, const fq_zech_mpoly_ctx_t ctx);

slong fq_zech_mpoly_get_term_var_exp_si(const fq_zech_mpoly_t A,
                            slong i, slong var, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_term_exp_fmpz(fq_zech_mpoly_t A,
                   slong i, fmpz * const * exp, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_set_term_exp_ui(fq_zech_mpoly_t A,
                    slong i, const ulong * exp, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term(fq_zech_mpoly_t M, const fq_zech_mpoly_t A,
                                       slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term_monomial(fq_zech_mpoly_t M,
              const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_push_term_fq_zech_fmpz(fq_zech_mpoly_t A,
         const fq_zech_t c, fmpz * const * exp, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_push_term_fq_zech_ui(fq_zech_mpoly_t A,
          const fq_zech_t c, const ulong * exp, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_sort_terms(fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_combine_like_terms(fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_reverse(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_assert_canonical(const fq_zech_mpoly_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_radix_sort1(fq_zech_mpoly_t A, slong left,
                 slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

void _fq_zech_mpoly_radix_sort(fq_zech_mpoly_t A, slong left,
                       slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask);

void _fq_zech_mpoly_push_exp_ffmpz(fq_zech_mpoly_t A,
                              const fmpz * exp, const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_push_exp_pfmpz(fq_zech_mpoly_t A,
                            fmpz * const * exp, const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_push_exp_ui(fq_zech_mpoly_t A,
                             const ulong * exp, const fq_zech_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

void fq_zech_mpoly_randtest_bound(fq_zech_mpoly_t A, flint_rand_t state,
                 slong length, ulong exp_bound, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_randtest_bounds(fq_zech_mpoly_t A, flint_rand_t state,
              slong length, ulong * exp_bounds, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_randtest_bits(fq_zech_mpoly_t A, flint_rand_t state,
            slong length, flint_bitcnt_t exp_bits, const fq_zech_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

slong _fq_zech_mpoly_add(
                         fq_zech_struct * coeff1,       ulong * exp1,
                         fq_zech_struct * coeff2, const ulong * exp2, slong len2,
                         fq_zech_struct * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx);

void fq_zech_mpoly_add_fq_zech(fq_zech_mpoly_t A,
                            const fq_zech_mpoly_t B, const fq_zech_t C,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_sub_fq_zech(fq_zech_mpoly_t A,
                            const fq_zech_mpoly_t B, const fq_zech_t C,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_add(fq_zech_mpoly_t A,
                            const fq_zech_mpoly_t B, const fq_zech_mpoly_t C,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_sub(fq_zech_mpoly_t A,
                            const fq_zech_mpoly_t B, const fq_zech_mpoly_t C,
                                                const fq_zech_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

void fq_zech_mpoly_neg(fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_scalar_mul_fq_zech(fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B, const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_make_monic(fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_scalar_addmul_fq_zech(fq_zech_mpoly_t A,
         const fq_zech_mpoly_t B, const fq_zech_mpoly_t C, const fq_zech_t d,
                                                const fq_zech_mpoly_ctx_t ctx);

/* Differentiation **********************************************************/

void fq_zech_mpoly_derivative(fq_zech_mpoly_t A,
            const fq_zech_mpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

void fq_zech_mpoly_evaluate_one_fq_zech(fq_zech_mpoly_t A,
                    const fq_zech_mpoly_t B, slong var, const fq_zech_t val,
                                                const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_eval_all_fq_zech(fq_zech_t eval,
               const fq_zech_struct * Acoeffs, const ulong * Aexps, slong Alen,
                    flint_bitcnt_t Abits, fq_zech_struct * const * alphas,
                            const mpoly_ctx_t mctx, const fq_zech_ctx_t fqctx);

void fq_zech_mpoly_evaluate_all_fq_zech(fq_zech_t ev, const fq_zech_mpoly_t A,
                 fq_zech_struct * const * vals, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_compose_fq_zech_poly(fq_zech_poly_t A,
                     const fq_zech_mpoly_t B, fq_zech_poly_struct * const * C,
                                                const fq_zech_mpoly_ctx_t ctx);

void _fq_zech_mpoly_compose_mat(fq_zech_mpoly_t A,
                            const fq_zech_mpoly_t B, const fmpz_mat_t M,
              const fq_zech_mpoly_ctx_t ctxB, const fq_zech_mpoly_ctx_t ctxAC);

int fq_zech_mpoly_compose_fq_zech_mpoly_geobucket(fq_zech_mpoly_t A,
              const fq_zech_mpoly_t B, fq_zech_mpoly_struct * const * C,
              const fq_zech_mpoly_ctx_t ctxB, const fq_zech_mpoly_ctx_t ctxAC);

int fq_zech_mpoly_compose_fq_zech_mpoly_horner(fq_zech_mpoly_t A,
              const fq_zech_mpoly_t B, fq_zech_mpoly_struct * const * C,
              const fq_zech_mpoly_ctx_t ctxB, const fq_zech_mpoly_ctx_t ctxAC);

int fq_zech_mpoly_compose_fq_zech_mpoly(fq_zech_mpoly_t A,
              const fq_zech_mpoly_t B, fq_zech_mpoly_struct * const * C,
              const fq_zech_mpoly_ctx_t ctxB, const fq_zech_mpoly_ctx_t ctxAC);

void fq_zech_mpoly_compose_fq_zech_mpoly_gen(fq_zech_mpoly_t A,
                             const fq_zech_mpoly_t B, const slong * c,
              const fq_zech_mpoly_ctx_t ctxB, const fq_zech_mpoly_ctx_t ctxAC);


/* Multiplication ************************************************************/

void fq_zech_mpoly_mul(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                       const fq_zech_mpoly_t C, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_mul_johnson(fq_zech_mpoly_t poly1,
                    const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                const fq_zech_mpoly_ctx_t ctx);

slong _fq_zech_mpoly_mul_johnson(
                    fq_zech_struct ** coeff1, ulong ** exp1, slong * alloc,
             const fq_zech_struct * coeff2, const ulong * exp2, slong len2,
             const fq_zech_struct * coeff3, const ulong * exp3, slong len3,
  flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx);


/* Powering ******************************************************************/

int fq_zech_mpoly_pow_fmpz(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                const fmpz_t k, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_pow_ui(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                       ulong k, const fq_zech_mpoly_ctx_t ctx);


/* Division ******************************************************************/

int fq_zech_mpoly_divides(fq_zech_mpoly_t Q,
                         const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_div(fq_zech_mpoly_t Q,
                            const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_divrem(fq_zech_mpoly_t Q, fq_zech_mpoly_t R,
                            const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_divrem_ideal(fq_zech_mpoly_struct ** Q,
                                  fq_zech_mpoly_t R, const fq_zech_mpoly_t A,
                                  fq_zech_mpoly_struct * const * B, slong len,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_divides_monagan_pearce(fq_zech_mpoly_t poly1,
                  const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_div_monagan_pearce(fq_zech_mpoly_t q,
                      const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_divrem_monagan_pearce(fq_zech_mpoly_t q, fq_zech_mpoly_t r,
                      const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_divrem_ideal_monagan_pearce(
                        fq_zech_mpoly_struct ** q, fq_zech_mpoly_t r,
            const fq_zech_mpoly_t poly2, fq_zech_mpoly_struct * const * poly3,
                                      slong len, const fq_zech_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

int fq_zech_mpoly_gcd(fq_zech_mpoly_t G, const fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx);

int _fq_zech_mpoly_gcd(fq_zech_mpoly_t G, flint_bitcnt_t Gbits,
                            const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_gcd_cofactors(fq_zech_mpoly_t G,
          fq_zech_mpoly_t Abar, fq_zech_mpoly_t Bbar, const fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx);

int _fq_zech_mpoly_gcd_cofactors(
                                fq_zech_mpoly_t G, flint_bitcnt_t Gbits,
                                fq_zech_mpoly_t Abar, flint_bitcnt_t Abarbits,
                                fq_zech_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
                             const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_gcd_brown(fq_zech_mpoly_t G,
                            const fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_gcd_zippel(fq_zech_mpoly_t G, const fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_deflation(fmpz * shift, fmpz * stride,
                       const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_deflate(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_inflate(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
       const fmpz * shift, const fmpz * stride, const fq_zech_mpoly_ctx_t ctx);


/* Univariates ***************************************************************/

void fq_zech_mpoly_univar_init(fq_zech_mpoly_univar_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_univar_clear(fq_zech_mpoly_univar_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_univar_fit_length(fq_zech_mpoly_univar_t A,
                                  slong length, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_univar_print_pretty(const fq_zech_mpoly_univar_t A,
                               const char ** x, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_univar_assert_canonical(fq_zech_mpoly_univar_t A,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_to_univar(fq_zech_mpoly_univar_t A,
            const fq_zech_mpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_from_univar_bits(fq_zech_mpoly_t A, flint_bitcnt_t Abits,
     const fq_zech_mpoly_univar_t B, slong var, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_from_univar(fq_zech_mpoly_t A,
     const fq_zech_mpoly_univar_t B, slong var, const fq_zech_mpoly_ctx_t ctx);

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_univar_swap(fq_zech_mpoly_univar_t A,
                       fq_zech_mpoly_univar_t B, const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fq_zech_mpoly_univar_struct, *A, *B);
}

FQ_ZECH_MPOLY_INLINE
int fq_zech_mpoly_univar_degree_fits_si(const fq_zech_mpoly_univar_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_univar_length(const fq_zech_mpoly_univar_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    return A->length;
}

FQ_ZECH_MPOLY_INLINE
slong fq_zech_mpoly_univar_get_term_exp_si(fq_zech_mpoly_univar_t A, slong i,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_univar_get_term_coeff(fq_zech_mpoly_t c,
        const fq_zech_mpoly_univar_t A, slong i, const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fq_zech_mpoly_set(c, A->coeffs + i, ctx);
}

FQ_ZECH_MPOLY_INLINE
void fq_zech_mpoly_univar_swap_term_coeff(fq_zech_mpoly_t c,
              fq_zech_mpoly_univar_t A, slong i, const fq_zech_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fq_zech_mpoly_swap(c, A->coeffs + i, ctx);
}


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/


int _fq_zech_mpoly_get_nmod_mpoly(
    nmod_mpoly_t s,
    const nmod_mpoly_ctx_t sctx,
    const fq_zech_mpoly_t t,
    const fq_zech_mpoly_ctx_t tctx);

void _fq_zech_mpoly_set_nmod_mpoly(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx);

void fq_zech_mpolyl_lead_coeff(
    fq_zech_mpoly_t c,
    const fq_zech_mpoly_t A,
    slong num_vars,
    const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_pow_rmul(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                       ulong k, const fq_zech_mpoly_ctx_t ctx);

int fq_zech_mpoly_repack_bits_inplace(fq_zech_mpoly_t A,
                         flint_bitcnt_t Abits, const fq_zech_mpoly_ctx_t ctx);


void fq_zech_mpoly_ctx_change_modulus(fq_zech_mpoly_ctx_t ctx,
                                                                    slong deg);

void _fq_zech_mpoly_get_fq_nmod_mpoly(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctxA,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_ctx_t ctxB);

void _fq_zech_mpoly_set_fq_nmod_mpoly(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctxA,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctxB);


/* mpolyu ********************************************************************/

int fq_zech_mpolyu_is_canonical(const fq_zech_mpolyu_t poly,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyu_init(fq_zech_mpolyu_t A, flint_bitcnt_t bits,
                                                const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyu_clear(fq_zech_mpolyu_t A,
                                               const fq_zech_mpoly_ctx_t uctx);

void fq_zech_mpolyu_swap(fq_zech_mpolyu_t A, fq_zech_mpolyu_t B);

void fq_zech_mpolyu_zero(fq_zech_mpolyu_t A,
                                               const fq_zech_mpoly_ctx_t uctx);

int fq_zech_mpolyu_is_one(fq_zech_mpolyu_t A,
                                               const fq_zech_mpoly_ctx_t uctx);

void fq_zech_mpolyu_print_pretty(const fq_zech_mpolyu_t poly,
                               const char ** x, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpolyu_fit_length(fq_zech_mpolyu_t A, slong length,
                                               const fq_zech_mpoly_ctx_t uctx);

void fq_zech_mpolyu_one(fq_zech_mpolyu_t A,
                                               const fq_zech_mpoly_ctx_t uctx);

#ifdef __cplusplus
}
#endif

#endif

