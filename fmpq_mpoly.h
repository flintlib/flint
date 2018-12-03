/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MPOLY_H
#define FMPQ_MPOLY_H

#ifdef FMPQ_MPOLY_INLINES_C
#define FMPQ_MPOLY_INLINE FLINT_DLL
#else
#define FMPQ_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "fmpq_poly.h"
#include "fmpz_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Context object ************************************************************/

typedef struct
{
    fmpz_mpoly_ctx_t zctx;
} fmpq_mpoly_ctx_struct;

typedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1];

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

/*
    A polynomial f is represented as
        content * zpoly,
    where zpoly should have positive leading coefficient and trivial content.
    If f is zero, then the representation should have
        content = 0 and zpoly = 0
*/

typedef struct
{                       /* non zero case:                   |  zero case: */
    fmpq_t content;     /* positive or negative content     |  (or zero)  */
    fmpz_mpoly_t zpoly; /* contentless poly, lc is positive |  (or zero)  */
} fmpq_mpoly_struct;

typedef fmpq_mpoly_struct fmpq_mpoly_t[1];


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
void fmpq_mpoly_init3(fmpq_mpoly_t A, slong alloc, mp_bitcnt_t bits,
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
                                  mp_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_bits(A->zpoly, bits, ctx->zctx);
}


/* Input/output **************************************************************/

FLINT_DLL int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char * str,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_fprint_pretty(FILE * file, 
            const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_print_pretty(const fmpq_mpoly_t A,
                                   const char ** x, const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_mpoly_fprint_pretty(stdout, A, x, ctx);
}


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
    fmpq_mpoly_struct t = *A;
    *A = *B;
    *B = t;
}


/* Constants *****************************************************************/

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_is_fmpz(A->zpoly, ctx->zctx);
}

FLINT_DLL void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_fmpq(fmpq_mpoly_t A,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_fmpz(fmpq_mpoly_t A,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_ui(fmpq_mpoly_t A,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_si(fmpq_mpoly_t A,
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

FLINT_DLL int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t A, const fmpq_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t A, const fmpz_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_ui(const fmpq_mpoly_t A,   ulong        c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_si(const fmpq_mpoly_t A,   slong        c,
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


/* Coefficients **************************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_denominator(fmpz_t d, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set(d, fmpq_denref(A->content));
}

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c,
                              const fmpq_mpoly_t A, const fmpq_mpoly_t M,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t A,
                                    const fmpq_t c, const fmpq_mpoly_t M,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A,
                 const fmpq_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A,
                const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t A,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A,
                                 const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A,
                               fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t A,
                                const ulong * exp, const fmpq_mpoly_ctx_t ctx);


/* container operations ******************************************************/

FLINT_DLL int fmpq_mpoly_is_canonical(const fmpq_mpoly_t A,
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

FLINT_DLL void fmpq_mpoly_get_term_coeff_fmpq(fmpq_t c, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_term_coeff_fmpq(fmpq_mpoly_t A,
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

FMPQ_MPOLY_INLINE
void fmpq_mpoly_get_term_exp_fmpz(fmpz ** exps, const fmpq_mpoly_t A,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_get_term_exp_fmpz(exps, A->zpoly, i, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_get_term_exp_ui(ulong * exps, const fmpq_mpoly_t A,
                                          slong i, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_get_term_exp_ui(exps, A->zpoly, i, ctx->zctx);
}

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

FLINT_DLL void fmpq_mpoly_push_term_fmpq_fmpz(fmpq_mpoly_t A,
               const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_fmpz_fmpz(fmpq_mpoly_t A,
               const fmpz_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_ui_fmpz(fmpq_mpoly_t A,
                      ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_si_fmpz(fmpq_mpoly_t A,
                      slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_t A,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_fmpz_ui(fmpq_mpoly_t A,
                const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_ui_ui(fmpq_mpoly_t A,
                       ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_push_term_si_ui(fmpq_mpoly_t A,
                       slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_reduce(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx);

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

FLINT_DLL void fmpq_mpoly_reverse(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_emplacebackterm_fmpq_ui(fmpq_mpoly_t poly,
                       fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_emplacebackterm_fmpq_fmpz(fmpq_mpoly_t poly,
                     fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bounds(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, ulong * exp_bounds,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bounds(A->zpoly, state,
                                     length, coeff_bits, exp_bounds, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bound(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, ulong exp_bound,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bound(A->zpoly, state,
                                     length, coeff_bits, exp_bound, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bits(fmpq_mpoly_t A, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, mp_bitcnt_t exp_bits,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bits(A->zpoly, state,
                                      length, coeff_bits, exp_bits, ctx->zctx);
    fmpq_randtest_not_zero(A->content, state, coeff_bits + 1);
    fmpq_mpoly_reduce(A, ctx);
}


/* Basic arithmetic **********************************************************/

FLINT_DLL void fmpq_mpoly_add_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          slong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          slong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_neg(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_neg(A->content, B->content);
    fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
}

FLINT_DLL void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t A,
             const fmpq_mpoly_t B, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_make_monic(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);


/* Differentiation/Integration ***********************************************/

FLINT_DLL void fmpq_mpoly_derivative(fmpq_mpoly_t A,
                  const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_integral(fmpq_mpoly_t A,
                  const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

FLINT_DLL void _fmpq_mpoly_rescale(fmpq_t Acontent, fmpz * Acoeff,
       const fmpq_mpoly_t B, const fmpq * scales,  const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A,
                              fmpq * const * vals, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A,
                           const fmpq_mpoly_t B, slong var, const fmpq_t val,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A,
                         const fmpq_mpoly_t B, fmpq_poly_struct * const * C,
                                                  const fmpq_mpoly_ctx_t ctxB);

FLINT_DLL void fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A,
                   const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C,
                    const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC);



/* Multiplication ************************************************************/

FLINT_DLL void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                             const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL void fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                   const fmpz_t k, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                          ulong k, const fmpq_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL int fmpq_mpoly_divides(fmpq_mpoly_t poly1,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_div(fmpq_mpoly_t q,
                     const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_divrem(fmpq_mpoly_t q, fmpq_mpoly_t r,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** q, fmpq_mpoly_t r,
    const fmpq_mpoly_t poly2, fmpq_mpoly_struct * const * poly3, slong len,
                                                   const fmpq_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

FLINT_DLL void fmpq_mpoly_term_content(fmpq_mpoly_t M, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A,
                             const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx);
/*
not implemented yet
FLINT_DLL void fmpq_mpoly_resultant(fmpq_mpoly_t poly1,
                const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_discriminant(fmpq_mpoly_t poly1,
              const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx);
*/

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

/* geobuckets ****************************************************************/
typedef struct fmpq_mpoly_geobucket
{
    fmpq_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fmpq_mpoly_geobucket_struct;

typedef fmpq_mpoly_geobucket_struct fmpq_mpoly_geobucket_t[1];

FLINT_DLL void fmpq_mpoly_geobucket_init(fmpq_mpoly_geobucket_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_clear(fmpq_mpoly_geobucket_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_empty(fmpq_mpoly_t p,
                         fmpq_mpoly_geobucket_t B, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_print(fmpq_mpoly_geobucket_t B,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_fit_length(fmpq_mpoly_geobucket_t B,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_geobucket_fix(fmpq_mpoly_geobucket_t B, slong i,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_add(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_sub(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set_fmpz(fmpq_mpoly_geobucket_t B,
                                         fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set_fmpq(fmpq_mpoly_geobucket_t B,
                                         fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_gen(fmpq_mpoly_geobucket_t B, slong var,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_add_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_sub_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_neg_inplace(fmpq_mpoly_geobucket_t B1,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_mul_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_pow_ui_inplace(fmpq_mpoly_geobucket_t B1,
                                          ulong k, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_pow_fmpz_inplace(fmpq_mpoly_geobucket_t B1,
                                   const fmpz_t k, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_geobucket_divides_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
FMPQ_MPOLY_INLINE
void fmpq_mpoly_remainder_test(const fmpq_mpoly_t r, const fmpq_mpoly_t g,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_remainder_strongtest(r->zpoly, g->zpoly, ctx->zctx);
}



#ifdef __cplusplus
}
#endif

#endif
