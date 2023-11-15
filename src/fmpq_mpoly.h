/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MPOLY_H
#define FMPQ_MPOLY_H

#ifdef FMPQ_MPOLY_INLINES_C
#define FMPQ_MPOLY_INLINE
#else
#define FMPQ_MPOLY_INLINE static inline
#endif

#include "fmpq.h"
#include "fmpq_types.h"
#include "fmpz_mpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Context object ************************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx,
                                            slong nvars, const ordering_t ord)
{
    fmpz_mpoly_ctx_init(ctx->zctx, nvars, ord);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_init_rand(fmpq_mpoly_ctx_t ctx,
                                           flint_rand_t state, slong max_nvars)
{
    fmpz_mpoly_ctx_init_rand(ctx->zctx, state, max_nvars);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_ctx_clear(ctx->zctx);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_ctx_nvars(const fmpq_mpoly_ctx_t ctx)
{
    return ctx->zctx->minfo->nvars;
}

FMPQ_MPOLY_INLINE
ordering_t fmpq_mpoly_ctx_ord(const fmpq_mpoly_ctx_t ctx)
{
    return ctx->zctx->minfo->ord;
}

/* Polynomials over Q ********************************************************/

FMPQ_MPOLY_INLINE
fmpq * fmpq_mpoly_content_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->content;
}

FMPQ_MPOLY_INLINE
fmpz_mpoly_struct * fmpq_mpoly_zpoly_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly;
}

FMPQ_MPOLY_INLINE
fmpz * fmpq_mpoly_zpoly_term_coeff_ref(fmpq_mpoly_t A, slong i,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < A->zpoly->length);
    return A->zpoly->coeffs + i;
}


/* Internal type definitions *************************************************/

/*
    fmpq_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fmpq_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fmpq_mpoly_univar_struct;

typedef fmpq_mpoly_univar_struct fmpq_mpoly_univar_t[1];


/*  Memory management ********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_init(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(A->content);
    fmpz_mpoly_init(A->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_init2(fmpq_mpoly_t A, slong alloc,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(A->content);
    fmpz_mpoly_init2(A->zpoly, alloc, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_init3(fmpq_mpoly_t A, slong alloc, flint_bitcnt_t bits,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(A->content);
    fmpz_mpoly_init3(A->zpoly, alloc, bits, ctx->zctx);
}


FMPQ_MPOLY_INLINE
void fmpq_mpoly_realloc(fmpq_mpoly_t A, slong alloc,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_realloc(A->zpoly, alloc, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_fit_length(fmpq_mpoly_t A, slong len,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_length(A->zpoly, len, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_clear(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_clear(A->zpoly, ctx->zctx);
    fmpq_clear(A->content);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_fit_bits(fmpq_mpoly_t A,
                                  flint_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_bits(A->zpoly, bits, ctx->zctx);
}


/* Input/output **************************************************************/

int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char * str,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fmpq_mpoly_fprint_pretty(FILE * file, const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx);
#endif

int fmpq_mpoly_print_pretty(const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx);

/*  Basic manipulation *******************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_gen(fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_one(A->content);
    fmpz_mpoly_gen(A->zpoly, var, ctx->zctx);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_gen(const fmpq_mpoly_t A,
                                         slong var, const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_is_one(A->content)
          && fmpz_mpoly_is_gen(A->zpoly, var, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_set(A->content, B->content);
    fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_equal(const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_equal(A->content, B->content)
        && fmpz_mpoly_equal(A->zpoly, B->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_swap(fmpq_mpoly_t A,
                                fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fmpq_mpoly_struct, *A, *B);
}


/* Constants *****************************************************************/

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_is_fmpz(A->zpoly, ctx->zctx);
}

void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_fmpq(fmpq_mpoly_t A,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_fmpz(fmpq_mpoly_t A,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_ui(fmpq_mpoly_t A,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_si(fmpq_mpoly_t A,
                                          slong c, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
void fmpq_mpoly_zero(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_zero(A->content);
    fmpz_mpoly_zero(A->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_one(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_one(A->content);
    fmpz_mpoly_one(A->zpoly, ctx->zctx);
}

int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t A, const fmpq_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t A, const fmpz_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_equal_ui(const fmpq_mpoly_t A,   ulong        c,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_equal_si(const fmpq_mpoly_t A,   slong        c,
                                                   const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_zero(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->length == WORD(0);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_one(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
   return fmpq_is_one(A->content)
          && fmpz_mpoly_is_one(A->zpoly, ctx->zctx);
}


/* Degrees *******************************************************************/

FMPQ_MPOLY_INLINE
int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_degrees_fit_si(A->zpoly->exps,
                                       A->zpoly->length, A->zpoly->bits,
                                                             ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t A, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->zpoly->exps, A->zpoly->length,
                                        A->zpoly->bits, var, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_degree_si(const fmpq_mpoly_t A, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->zpoly->exps, A->zpoly->length,
                                        A->zpoly->bits, var, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_total_degree_fits_si(const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(tdeg, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_total_degree_si(const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_used_vars(int * used, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_used_vars(used, A->zpoly, ctx->zctx);
}

/* Coefficients **************************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_get_denominator(fmpz_t d, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set(d, fmpq_denref(A->content));
}

int fmpq_mpoly_is_monic(const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c,
                              const fmpq_mpoly_t A, const fmpq_mpoly_t M,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t A,
                                    const fmpq_t c, const fmpq_mpoly_t M,
                                                   const fmpq_mpoly_ctx_t ctx);

void _fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A,
                 const fmpq_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A,
                const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t A,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void _fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A,
                                 const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A,
                               fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t A,
                                const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_coeff_vars_ui(fmpq_mpoly_t C,
                 const fmpq_mpoly_t A, const slong * vars, const ulong * exps,
                                     slong length, const fmpq_mpoly_ctx_t ctx);

/* conversion ****************************************************************/

int fmpq_mpoly_is_fmpq_poly(const fmpq_mpoly_t A,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_get_fmpq_poly(fmpq_poly_t A,  const fmpq_mpoly_t B,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_fmpq_poly(fmpq_mpoly_t A, const fmpq_poly_t B,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

/* comparison ****************************************************************/

int fmpq_mpoly_cmp(const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);


/* container operations ******************************************************/

int fmpq_mpoly_is_canonical(const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_length(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->length;
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_resize(fmpq_mpoly_t A, slong new_length,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_resize(A->zpoly, new_length, ctx->zctx);
}

void fmpq_mpoly_get_term_coeff_fmpq(fmpq_t c, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_set_term_coeff_fmpq(fmpq_mpoly_t A,
                          slong i, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_term_exp_fits_ui(const fmpq_mpoly_t A,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_term_exp_fits_ui(A->zpoly->exps,
                                          A->zpoly->bits, i, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_term_exp_fits_si(const fmpq_mpoly_t A,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_term_exp_fits_si(A->zpoly->exps,
                                          A->zpoly->bits, i, ctx->zctx->minfo);
}

void fmpq_mpoly_get_term_exp_fmpz(fmpz ** exps, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_term_exp_ui(ulong * exps, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_term_exp_si(slong * exps, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

ulong fmpq_mpoly_get_term_var_exp_ui(const fmpq_mpoly_t A,
                               slong i, slong var, const fmpq_mpoly_ctx_t ctx);

slong fmpq_mpoly_get_term_var_exp_si(const fmpq_mpoly_t A,
                               slong i, slong var, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set_term_exp_fmpz(fmpq_mpoly_t A,
                      slong i, fmpz * const * exps, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_term_exp_fmpz(A->zpoly, i, exps, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set_term_exp_ui(fmpq_mpoly_t A,
                       slong i, const ulong * exps, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_term_exp_ui(A->zpoly, i, exps, ctx->zctx);
}

void fmpq_mpoly_get_term(fmpq_mpoly_t M, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_get_term_monomial(fmpq_mpoly_t M, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpq_fmpz(fmpq_mpoly_t A,
               const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpq_ffmpz(fmpq_mpoly_t A, const fmpq_t c,
                                  const fmpz *exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpz_fmpz(fmpq_mpoly_t A, const fmpz_t c,
                                        fmpz *const *exp,
                                        const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpz_ffmpz(fmpq_mpoly_t A, const fmpz_t c,
                                     const fmpz *exp,
                                     const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_ui_fmpz(fmpq_mpoly_t A, ulong c, fmpz *const *exp,
                                      const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_ui_ffmpz(fmpq_mpoly_t A, ulong c, const fmpz *exp,
                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_si_fmpz(fmpq_mpoly_t A, slong c, fmpz *const *exp,
                                  const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_si_ffmpz(fmpq_mpoly_t A, slong c, const fmpz *exp,
                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_t A,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_fmpz_ui(fmpq_mpoly_t A,
                const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_ui_ui(fmpq_mpoly_t A,
                       ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_push_term_si_ui(fmpq_mpoly_t A,
                       slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_reduce(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_reduce_easy(fmpq_mpoly_t A, slong easy_length,
                                                       const fmpq_mpoly_ctx_t);

FMPQ_MPOLY_INLINE
void fmpq_mpoly_sort_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_sort_terms(A->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_combine_like_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_combine_like_terms(A->zpoly, ctx->zctx);
    fmpq_mpoly_reduce(A, ctx);
}

void fmpq_mpoly_reverse(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx);

void _fmpq_mpoly_push_rescale(fmpq_mpoly_t A,
                                         fmpq_t C, const fmpq_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bounds(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, flint_bitcnt_t coeff_bits, ulong * exp_bounds,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bounds(A->zpoly, state,
                                     length, coeff_bits, exp_bounds, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bound(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, flint_bitcnt_t coeff_bits, ulong exp_bound,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bound(A->zpoly, state,
                                     length, coeff_bits, exp_bound, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bits(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bits(A->zpoly, state,
                                      length, coeff_bits, exp_bits, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}


/* Basic arithmetic **********************************************************/

void fmpq_mpoly_add_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_add_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_add_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_add_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          slong c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_sub_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_sub_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          slong c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_neg(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_neg(A->content, B->content);
    fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
}

void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, ulong        c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, slong        c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, ulong        c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, slong        c, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_make_monic(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

void _fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);


/* Differentiation/Integration ***********************************************/

void fmpq_mpoly_derivative(fmpq_mpoly_t A,
                  const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_integral(fmpq_mpoly_t A,
                  const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

int _fmpq_mpoly_rescale(fmpq_t Acontent, fmpz * Acoeff,
       const fmpq_mpoly_t B, const fmpq * scales,  const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A,
                              fmpq * const * vals, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A,
                           const fmpq_mpoly_t B, slong var, const fmpq_t val,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A,
                         const fmpq_mpoly_t B, fmpq_poly_struct * const * C,
                                                  const fmpq_mpoly_ctx_t ctxB);

int fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A,
                   const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C,
                    const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC);

void fmpq_mpoly_compose_fmpq_mpoly_gen(fmpq_mpoly_t A,
                     const fmpq_mpoly_t B, const slong * c,
                    const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC);


/* Multiplication ************************************************************/

void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

int fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t k, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong k, const fmpq_mpoly_ctx_t ctx);

/* Division ******************************************************************/

int fmpq_mpoly_divides(fmpq_mpoly_t poly1,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_div(fmpq_mpoly_t q,
                     const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_divrem(fmpq_mpoly_t q, fmpq_mpoly_t r,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** q, fmpq_mpoly_t r,
    const fmpq_mpoly_t poly2, fmpq_mpoly_struct * const * poly3, slong len,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Square root ***************************************************************/

int fmpq_mpoly_sqrt(fmpq_mpoly_t Q, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_square(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return fmpz_is_square(fmpq_numref(A->content)) &&
           fmpz_is_square(fmpq_denref(A->content)) &&
           fmpz_mpoly_is_square(A->zpoly, ctx->zctx);
}

/* GCD ***********************************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_content(fmpq_t g, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpq_abs(g, A->content);
}

void fmpq_mpoly_term_content(fmpq_mpoly_t M, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_content_vars(fmpq_mpoly_t g, const fmpq_mpoly_t A,
                  slong * vars, slong vars_length, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A,
                             const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_inflate(fmpq_mpoly_t A, const fmpq_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_cofactors(fmpq_mpoly_t G, fmpq_mpoly_t Abar,
             fmpq_mpoly_t Bbar, const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_hensel(fmpq_mpoly_t G,
       const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_brown(fmpq_mpoly_t G,
       const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_subresultant(fmpq_mpoly_t G,
       const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_zippel(fmpq_mpoly_t G,
       const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_gcd_zippel2(fmpq_mpoly_t G,
       const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_resultant(fmpq_mpoly_t R, const fmpq_mpoly_t A,
                  const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_discriminant(fmpq_mpoly_t R, const fmpq_mpoly_t A,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

void mpoly_void_ring_init_fmpq_mpoly_ctx(mpoly_void_ring_t R,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_repack_bits(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             flint_bitcnt_t Abits, const fmpq_mpoly_ctx_t ctx);

/* Univariates ***************************************************************/

void fmpq_mpoly_univar_init(fmpq_mpoly_univar_t A,
                                                const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_univar_clear(fmpq_mpoly_univar_t A,
                                                const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_univar_fit_length(fmpq_mpoly_univar_t A,
                                  slong length, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_univar_print_pretty(const fmpq_mpoly_univar_t A,
                               const char ** x, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_univar_assert_canonical(fmpq_mpoly_univar_t A,
                                                const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_to_univar(fmpq_mpoly_univar_t A,
            const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_from_univar_bits(fmpq_mpoly_t A, flint_bitcnt_t Abits,
     const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_from_univar(fmpq_mpoly_t A,
     const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
void fmpq_mpoly_univar_swap(fmpq_mpoly_univar_t A, fmpq_mpoly_univar_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fmpq_mpoly_univar_struct, *A, *B);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_univar_degree_fits_si(const fmpq_mpoly_univar_t A,
                                                 const fmpq_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_univar_length(const fmpq_mpoly_univar_t A,
                                                 const fmpq_mpoly_ctx_t ctx)
{
    return A->length;
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_univar_get_term_exp_si(fmpq_mpoly_univar_t A, slong i,
                                                 const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_univar_get_term_coeff(fmpq_mpoly_t c,
        const fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpq_mpoly_set(c, A->coeffs + i, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_univar_swap_term_coeff(fmpq_mpoly_t c,
              fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpq_mpoly_swap(c, A->coeffs + i, ctx);
}


/******************************************************************************

   Internal consistency checks

******************************************************************************/

void fmpq_mpoly_remainder_test(const fmpq_mpoly_t r, const fmpq_mpoly_t g, const fmpq_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
