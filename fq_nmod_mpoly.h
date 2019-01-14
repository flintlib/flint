/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_H
#define FQ_NMOD_MPOLY_H

#ifdef FQ_NMOD_MPOLY_INLINES_C
#define FQ_NMOD_MPOLY_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"

#include "nmod_mpoly.h"


#ifdef __cplusplus
 extern "C" {
#endif

/* Context object ************************************************************/

FLINT_DLL void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                                       mp_limb_t p, slong deg);

FLINT_DLL void fq_nmod_mpoly_ctx_init_rand(fq_nmod_mpoly_ctx_t ctx,
                                       flint_rand_t state, slong max_nvars,
                                          mp_bitcnt_t p_bits, slong deg_bound);

FLINT_DLL void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx);

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

FLINT_DLL void fq_nmod_mpoly_init(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc,
                              mp_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_fit_length(fq_nmod_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);


FQ_NMOD_MPOLY_INLINE
void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen, 
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->alloc);
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

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_fit_bits(fq_nmod_mpoly_t A, slong bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    ulong * t;

    if (A->bits < bits)
    {
        if (A->alloc != 0)
        {
            N = mpoly_words_per_exp(bits, ctx->minfo);
            t = flint_malloc(N*A->alloc*sizeof(ulong));
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length,
                                                                   ctx->minfo);
            flint_free(A->exps);
            A->exps = t;
        }

        A->bits = bits;
    }
}


/* Input/output **************************************************************/

FLINT_DLL int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_fprint_pretty(FILE * file, 
      const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
   return fq_nmod_mpoly_fprint_pretty(stdout, A, x, ctx);
}


/*  Basic manipulation *******************************************************/

FLINT_DLL void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_equal(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpoly_struct t = *A;
   *A = *B;
   *B = t;
}


/* Constants *****************************************************************/

FLINT_DLL int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   _fq_nmod_mpoly_set_length(A, 0, ctx);
}

FLINT_DLL void fq_nmod_mpoly_one(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

FLINT_DLL int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A,
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
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
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


/* Coefficients **************************************************************/

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c,
                          const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A,
                                  const fq_nmod_t c, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c,
                             const fq_nmod_mpoly_t A, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c,
                                const fq_nmod_mpoly_t A, const ulong * exp,
                                               const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                       const fq_nmod_t c, const fmpz * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                     const fq_nmod_t c, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A,
                                      const fq_nmod_t c, const ulong * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C,
           const fq_nmod_mpoly_t A, slong * vars, ulong * exps, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* comparison ****************************************************************/

FLINT_DLL int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

FLINT_DLL int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

FLINT_DLL void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, slong new_length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_coeff_fq_nmod(fq_nmod_mpoly_t A, 
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

FLINT_DLL void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_exp_ui(ulong * exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A,
                            slong i, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A,
                   slong i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A,
                    slong i, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A,
         const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A,
          const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_assert_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_radix_sort1(fq_nmod_mpoly_t A, slong left,
                 slong right, mp_bitcnt_t pos, ulong cmpmask, ulong totalmask);

FLINT_DLL void _fq_nmod_mpoly_radix_sort(fq_nmod_mpoly_t A, slong left,
                       slong right, mp_bitcnt_t pos, slong N, ulong * cmpmask);

FLINT_DLL void _fq_nmod_mpoly_push_exp_ffmpz(fq_nmod_mpoly_t A,
                              const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_push_exp_pfmpz(fq_nmod_mpoly_t A,
                            fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_push_exp_ui(fq_nmod_mpoly_t A,
                             const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FLINT_DLL void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state,
                 slong length, ulong exp_bound, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state,
              slong length, ulong * exp_bounds, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state,
            slong length, mp_bitcnt_t exp_bits, const fq_nmod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

FLINT_DLL void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_add(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

FLINT_DLL void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);


/* Differentiation **********************************************************/

FLINT_DLL void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);


/* Multiplication ************************************************************/

FLINT_DLL void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                       const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t poly1,
                    const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_mul_johnson(
                    fq_nmod_struct ** coeff1, ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  mp_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);


/* Powering ******************************************************************/

FLINT_DLL void fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_pow_rmul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);


/* Division ******************************************************************/

FLINT_DLL int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q,
                         const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q,
                                  fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A,
                                  fq_nmod_mpoly_struct * const * B, slong len,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t poly1,
                  const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_div_monagan_pearce(fq_nmod_mpoly_t q,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_monagan_pearce(fq_nmod_mpoly_t q, fq_nmod_mpoly_t r,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_ideal_monagan_pearce(
                        fq_nmod_mpoly_struct ** q, fq_nmod_mpoly_t r,
            const fq_nmod_mpoly_t poly2, fq_nmod_mpoly_struct * const * poly3,
                                      slong len, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_divides_monagan_pearce(
                  fq_nmod_struct ** coeff1,      ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  mp_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/


FLINT_DLL int fq_nmod_mpoly_repack_bits(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                             mp_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx);

/* geobuckets ****************************************************************/
typedef struct fq_nmod_mpoly_geobucket
{
    fq_nmod_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fq_nmod_mpoly_geobucket_struct;

typedef fq_nmod_mpoly_geobucket_struct fq_nmod_mpoly_geobucket_t[1];

FLINT_DLL void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p,
                   fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_print(fq_nmod_mpoly_geobucket_t B,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set_ui(fq_nmod_mpoly_geobucket_t B,
                                       ulong c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set_fq_nmod_gen(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_gen(fq_nmod_mpoly_geobucket_t B, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_add_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_sub_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_neg_inplace(fq_nmod_mpoly_geobucket_t B1,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_mul_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_pow_ui_inplace(fq_nmod_mpoly_geobucket_t B1,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_pow_fmpz_inplace(fq_nmod_mpoly_geobucket_t B1,
                                const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_geobucket_divides_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);


/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/
FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_remainder_strongtest(const fq_nmod_mpoly_t r,
                        const fq_nmod_mpoly_t g, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N, bits;
    ulong mask = 0;
    ulong * rexp, * gexp;

    bits = FLINT_MAX(r->bits, g->bits);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0 )
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0 )
        return;

    rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
    gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
    mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
    mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

    /* mask with high bit set in each field of exponent vector */
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_printf("fq_nmod_mpoly_remainder_strongtest FAILED i = %wd\n", i);
            flint_printf("rem ");fq_nmod_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fq_nmod_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
    }

    flint_free(rexp);
    flint_free(gexp);
}


#ifdef __cplusplus
}
#endif

#endif
