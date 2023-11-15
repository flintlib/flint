/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_H
#define FMPZ_MOD_MPOLY_H

#ifdef FMPZ_MOD_MPOLY_INLINES_C
#define FMPZ_MOD_MPOLY_INLINE
#else
#define FMPZ_MOD_MPOLY_INLINE static inline
#endif

#include "nmod_types.h"
#include "fmpz_mod.h"
#include "mpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* type definitions **********************************************************/

/*
    fmpz_mod_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fmpz_mod_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fmpz_mod_mpoly_univar_struct;

typedef fmpz_mod_mpoly_univar_struct fmpz_mod_mpoly_univar_t[1];


/* Context object ************************************************************/

void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx,
                      slong nvars, const ordering_t ord, const fmpz_t modulus);

void fmpz_mod_mpoly_ctx_init_rand(fmpz_mod_mpoly_ctx_t ctx,
                    flint_rand_t state, slong max_nvars, const fmpz_t modulus);

void fmpz_mod_mpoly_ctx_init_rand_bits_prime(fmpz_mod_mpoly_ctx_t ctx,
                 flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits);

void fmpz_mod_mpoly_ctx_init_rand_bits(fmpz_mod_mpoly_ctx_t ctx,
                 flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits);

void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_ctx_nvars(const fmpz_mod_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FMPZ_MOD_MPOLY_INLINE
ordering_t fmpz_mod_mpoly_ctx_ord(const fmpz_mod_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}

FMPZ_MOD_MPOLY_INLINE
const fmpz * fmpz_mod_mpoly_ctx_modulus(const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_ctx_modulus(ctx->ffinfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_ctx_get_modulus(fmpz_t m, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_set(m, fmpz_mod_mpoly_ctx_modulus(ctx));
}

/*  Memory management ********************************************************/


FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_init(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_init2(fmpz_mod_mpoly_t A, slong alloc,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_init3(fmpz_mod_mpoly_t A, slong alloc,
                          flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_realloc(fmpz_mod_mpoly_t A,
                                  slong alloc, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_fit_length(fmpz_mod_mpoly_t A, slong length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_fit_length_fit_bits(fmpz_mod_mpoly_t A,
               slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_fit_length_reset_bits(fmpz_mod_mpoly_t A,
               slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void _fmpz_mod_mpoly_fit_length(
    fmpz ** coeffs,
    slong * coeffs_alloc,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length)
{
    if (length > *coeffs_alloc)
    {
        slong i, old_alloc = *coeffs_alloc;
        slong new_alloc = FLINT_MAX(length, old_alloc*2);
        *coeffs_alloc = new_alloc;
        *coeffs = (fmpz *) flint_realloc(*coeffs, new_alloc*sizeof(fmpz));
        for (i = old_alloc; i < new_alloc; i++)
            fmpz_init(*coeffs + i);
    }

    if (N*length > *exps_alloc)
    {
        *exps_alloc = FLINT_MAX(N*length, *exps_alloc*2);
        *exps = (ulong *) flint_realloc(*exps, *exps_alloc*sizeof(ulong));
    }
}

FMPZ_MOD_MPOLY_INLINE
void _fmpz_mod_mpoly_set_length(fmpz_mod_mpoly_t A, slong newlen,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, ctx->minfo)*newlen <= A->exps_alloc);

    A->length = newlen;
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_truncate(fmpz_mod_mpoly_t A, slong newlen,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }
}


/* Input/output **************************************************************/

int fmpz_mod_mpoly_set_str_pretty(fmpz_mod_mpoly_t A,
            const char * str, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fmpz_mod_mpoly_fprint_pretty(FILE * file, const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);
#endif

int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

/*  Basic manipulation *******************************************************/

void fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A,
                                    slong var, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_equal(const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t A, fmpz_mod_mpoly_t B,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fmpz_mod_mpoly_struct, *A, *B);
}

/* Constants *****************************************************************/

int fmpz_mod_mpoly_is_fmpz(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_fmpz_mod(fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_ui(fmpz_mod_mpoly_t A,
                                      ulong c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_t A,
                                      slong c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
   _fmpz_mod_mpoly_set_length(A, 0, ctx);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_one(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_set_si(A, 1, ctx);
}

int fmpz_mod_mpoly_equal_fmpz(const fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_equal_ui(const fmpz_mod_mpoly_t A,
                                      ulong c, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_equal_si(const fmpz_mod_mpoly_t A,
                                      slong c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_zero(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
   return A->length < 1;
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
   return fmpz_mod_mpoly_equal_si(A, 1, ctx);
}

/* Degrees *******************************************************************/

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_degrees_fit_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degrees_si(slong * degs, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mod_mpoly_t A, slong var,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_degree_si(const fmpz_mod_mpoly_t A, slong var,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_total_degree_fits_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_total_degree_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_used_vars(int * used, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}

/* Coefficients **************************************************************/

void fmpz_mod_mpoly_get_coeff_fmpz_monomial(fmpz_t c,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t M,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_fmpz_monomial(fmpz_mod_mpoly_t A,
                                    const fmpz_t c, const fmpz_mod_mpoly_t M,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_coeff_fmpz_fmpz(fmpz_t c,
                                const fmpz_mod_mpoly_t A, fmpz * const * exp,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_coeff_fmpz_ui(fmpz_t c,
                                const fmpz_mod_mpoly_t A, const ulong * exp,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A,
             const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A,
           const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_ui_fmpz(fmpz_mod_mpoly_t A,
                  ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_si_fmpz(fmpz_mod_mpoly_t A,
                  slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t A,
            const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_ui_ui(fmpz_mod_mpoly_t A,
                   ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_coeff_si_ui(fmpz_mod_mpoly_t A,
                   slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_coeff_vars_ui(fmpz_mod_mpoly_t C,
             const fmpz_mod_mpoly_t A, const slong * vars, const ulong * exps,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE fmpz * fmpz_mod_mpoly_leadcoeff(fmpz_mod_mpoly_t A)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* conversion ****************************************************************/

int fmpz_mod_mpoly_is_fmpz_mod_poly(const fmpz_mod_mpoly_t A,
                                    slong var, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_get_fmpz_mod_poly(fmpz_mod_poly_t A,
          const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_set_fmpz_mod_poly(fmpz_mod_mpoly_t A,
                    flint_bitcnt_t Abits, const fmpz * Bcoeffs, slong Blen,
                                    slong var, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_fmpz_mod_poly(fmpz_mod_mpoly_t A,
           const fmpz_mod_poly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);

/* comparison ****************************************************************/

int fmpz_mod_mpoly_cmp(const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

int fmpz_mod_mpoly_is_canonical(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length;
}

void fmpz_mod_mpoly_resize(fmpz_mod_mpoly_t A, slong new_length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A, slong i,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_term_coeff_ui(fmpz_mod_mpoly_t A, slong i,
                                      ulong c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_term_coeff_si(fmpz_mod_mpoly_t A, slong i,
                                      slong c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_term_exp_fits_ui(const fmpz_mod_mpoly_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_term_exp_fits_si(const fmpz_mod_mpoly_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_term_exp_si(slong * exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

ulong fmpz_mod_mpoly_get_term_var_exp_ui(const fmpz_mod_mpoly_t A,
                           slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx);

slong fmpz_mod_mpoly_get_term_var_exp_si(const fmpz_mod_mpoly_t A,
                           slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_term_exp_fmpz(fmpz_mod_mpoly_t A, slong i,
                           fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, slong i,
                            const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_term(fmpz_mod_mpoly_t M,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_get_term_monomial(fmpz_mod_mpoly_t M,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_fmpz_fmpz(fmpz_mod_mpoly_t A,
           const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_fmpz_ffmpz(fmpz_mod_mpoly_t A,
           const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_ui_fmpz(fmpz_mod_mpoly_t A,
                  ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_ui_ffmpz(fmpz_mod_mpoly_t A,
                  ulong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_si_fmpz(fmpz_mod_mpoly_t A,
                  slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_si_ffmpz(fmpz_mod_mpoly_t A,
                  slong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_fmpz_ui(fmpz_mod_mpoly_t A,
            const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_ui_ui(fmpz_mod_mpoly_t A,
                   ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_push_term_si_ui(fmpz_mod_mpoly_t A,
                   slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_combine_like_terms(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_reverse(fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_assert_canonical(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_radix_sort1(fmpz *, ulong *, slong left,
              slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

void _fmpz_mod_mpoly_radix_sort(fmpz *, ulong *, slong left,
                    slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask);

void _fmpz_mod_mpoly_push_exp_ffmpz(fmpz_mod_mpoly_t A,
                             const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_push_exp_pfmpz(fmpz_mod_mpoly_t A,
                           fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_push_exp_ui(fmpz_mod_mpoly_t A,
                            const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A,
                        flint_rand_t state, slong length, ulong * exp_bounds,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_t A,
                            flint_rand_t state, slong length, ulong exp_bound,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A,
                    flint_rand_t state, slong length, flint_bitcnt_t exp_bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

void fmpz_mod_mpoly_add_fmpz_mod(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_add_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_add_ui(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, ulong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_add_si(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, slong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_sub_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_sub_ui(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, ulong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_sub_si(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, slong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_add(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                     const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_sub(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                     const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx);

/* Scalar operations *********************************************************/

void fmpz_mod_mpoly_neg(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_scalar_mul_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_scalar_mul_ui(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, ulong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_scalar_mul_si(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, slong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_scalar_addmul_fmpz(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                               const fmpz_t d, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_make_monic(fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);


/* Differentiation ************************************************************/

void fmpz_mod_mpoly_derivative(fmpz_mod_mpoly_t A,
          const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

void _fmpz_mod_mpoly_eval_all_fmpz_mod(fmpz_t eval,
                        const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
                            flint_bitcnt_t Abits, const fmpz * alphas,
                            const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fctx);

void fmpz_mod_mpoly_evaluate_all_fmpz(fmpz_t eval,
                            const fmpz_mod_mpoly_t A, fmpz * const * alphas,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_evaluate_one_fmpz(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_t B, slong var, const fmpz_t val,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_compose_mat(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, const fmpz_mat_t M,
            const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC);

int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_geobucket(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C,
            const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC);

int fmpz_mod_mpoly_compose_fmpz_mod_mpoly(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C,
            const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC);

/* Multiplication ************************************************************/

void fmpz_mod_mpoly_mul(fmpz_mod_mpoly_t A,
                          const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_mul_johnson(fmpz_mod_mpoly_t A,
                          const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_mul_johnson_maxfields(fmpz_mod_mpoly_t A,
                                const fmpz_mod_mpoly_t B, fmpz * maxBfields,
                                const fmpz_mod_mpoly_t C, fmpz * maxCfields,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_mul_dense(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int _fmpz_mod_mpoly_mul_dense_maxfields(fmpz_mod_mpoly_t P,
                                const fmpz_mod_mpoly_t A, fmpz * maxAfields,
                                const fmpz_mod_mpoly_t B, fmpz * maxBfields,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

int fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_t A,
     const fmpz_mod_mpoly_t B, const fmpz_t k, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_pow_rmul(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx);

/* Division ******************************************************************/

int fmpz_mod_mpoly_divides(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_div(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_divrem(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_divrem_ideal(fmpz_mod_mpoly_struct ** Q,
                            fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A,
                            fmpz_mod_mpoly_struct * const * B, slong len,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_divexact(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A,
                      const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_divides(Q, A, B, ctx))
        return;

    flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_divexact: nonexact division");
}

int _fmpz_mod_mpoly_divides_dense_maxfields(fmpz_mod_mpoly_t Q,
                                const fmpz_mod_mpoly_t A, fmpz * maxAfields,
                                const fmpz_mod_mpoly_t B, fmpz * maxBfields,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_divides_dense(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int _fmpz_mod_mpoly_divides_monagan_pearce_maxfields(
            fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, fmpz * maxAfields,
                                const fmpz_mod_mpoly_t B, fmpz * maxBfields,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_divides_monagan_pearce(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_div_monagan_pearce(fmpz_mod_mpoly_t Q,
                          const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_divrem_monagan_pearce(fmpz_mod_mpoly_t Q,
        fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_divrem_ideal_monagan_pearce(
        fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R,
        const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Square root ***************************************************************/

int fmpz_mod_mpoly_sqrt_heap(fmpz_mod_mpoly_t Q,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_sqrt(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_sqrt_heap(Q, A, ctx);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_square(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    int res;
    fmpz_mod_mpoly_t Q;
    fmpz_mod_mpoly_init(Q, ctx);
    res = fmpz_mod_mpoly_sqrt_heap(Q, A, ctx);
    fmpz_mod_mpoly_clear(Q, ctx);
    return res;
}

int fmpz_mod_mpoly_quadratic_root(fmpz_mod_mpoly_t Q,
                    const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

void fmpz_mod_mpoly_term_content(fmpz_mod_mpoly_t M,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_content_vars(fmpz_mod_mpoly_t g,
                    const fmpz_mod_mpoly_t A, slong * vars, slong vars_length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_cofactors(fmpz_mod_mpoly_t G,
                        fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_subresultant(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_brown(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_hensel(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_zippel(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_gcd_zippel2(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_deflation(fmpz * shift, fmpz * stride,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_deflate(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_inflate(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Univariates ***************************************************************/

void fmpz_mod_mpoly_univar_init(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_univar_clear(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_univar_fit_length(fmpz_mod_mpoly_univar_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_univar_print_pretty(const fmpz_mod_mpoly_univar_t A,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_univar_assert_canonical(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_zero(fmpz_mod_mpoly_univar_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

void fmpz_mod_mpoly_univar_set_coeff_ui(fmpz_mod_mpoly_univar_t A,
            ulong e, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_to_univar(fmpz_mod_mpoly_univar_t A,
          const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A,
            flint_bitcnt_t Abits, const fmpz_mod_mpoly_univar_t B, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A,
                                  const fmpz_mod_mpoly_univar_t B, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_swap(fmpz_mod_mpoly_univar_t A,
                     fmpz_mod_mpoly_univar_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_SWAP(fmpz_mod_mpoly_univar_struct, *A, *B);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_univar_degree_fits_si(const fmpz_mod_mpoly_univar_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_univar_length(const fmpz_mod_mpoly_univar_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_univar_get_term_exp_si(fmpz_mod_mpoly_univar_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_get_term_coeff(fmpz_mod_mpoly_t c,
      const fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mod_mpoly_set(c, A->coeffs + i, ctx);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_swap_term_coeff(fmpz_mod_mpoly_t c,
            fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mod_mpoly_swap(c, A->coeffs + i, ctx);
}

int fmpz_mod_mpoly_univar_pseudo_gcd(fmpz_mod_mpoly_univar_t Gx,
        const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_univar_resultant(fmpz_mod_mpoly_t R,
        const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_univar_discriminant(fmpz_mod_mpoly_t D,
             const fmpz_mod_mpoly_univar_t Fx, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_resultant(fmpz_mod_mpoly_t R,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                    slong var, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_discriminant(fmpz_mod_mpoly_t R,
          const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

FMPZ_MOD_MPOLY_INLINE
void _fmpz_mod_mpoly_clear_dense_mock(fmpz_mod_poly_t D)
{
    flint_free(D->coeffs);
}

void _fmpz_mod_mpoly_init_dense_mock(fmpz_mod_poly_t D,
                        const fmpz_mod_mpoly_t A, const slong * Adeg_bounds,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void mpoly_void_ring_init_fmpz_mod_mpoly_ctx(mpoly_void_ring_t R,
                                              const fmpz_mod_mpoly_ctx_t ctx);

/* geobuckets ****************************************************************/

typedef struct fmpz_mod_mpoly_geobucket
{
    fmpz_mod_mpoly_struct polys[FLINT_BITS/2];
    fmpz_mod_mpoly_struct temps[FLINT_BITS/2];
    slong length;
} fmpz_mod_mpoly_geobucket_struct;

typedef fmpz_mod_mpoly_geobucket_struct fmpz_mod_mpoly_geobucket_t[1];

void fmpz_mod_mpoly_geobucket_init(fmpz_mod_mpoly_geobucket_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_clear(fmpz_mod_mpoly_geobucket_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_empty(fmpz_mod_mpoly_t p,
                 fmpz_mod_mpoly_geobucket_t B, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_fit_length(fmpz_mod_mpoly_geobucket_t B,
                                      slong i, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_set(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_add(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_geobucket_sub(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyl_lead_coeff(fmpz_mod_mpoly_t c,
     const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyl_content(fmpz_mod_mpoly_t g,
     const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(fmpz_mod_poly_t A,
                const fmpz_mod_mpoly_t B, slong var, const ulong * Bshift,
                        const ulong * Bstride, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(fmpz_mod_mpoly_t A,
                    flint_bitcnt_t Abits, const fmpz_mod_poly_t B, slong var,
                    const ulong * Ashift, const ulong * Astride,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_set_nmod_mpoly(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_ctx_t ctx, const nmod_mpoly_t nA,
                                                  const nmod_mpoly_ctx_t nctx);

void _fmpz_mod_mpoly_get_nmod_mpoly(nmod_mpoly_t nA,
                        const nmod_mpoly_ctx_t nctx, const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_repack_bits(fmpz_mod_mpoly_t A,
                            const fmpz_mod_mpoly_t B, flint_bitcnt_t Abits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_repack_bits_inplace(fmpz_mod_mpoly_t A,
                         flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_to_mpolyl_perm_deflate(fmpz_mod_mpoly_t A,
                    const fmpz_mod_mpoly_ctx_t lctx, const fmpz_mod_mpoly_t B,
                    const fmpz_mod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride);

void fmpz_mod_mpoly_from_mpolyl_perm_inflate(fmpz_mod_mpoly_t A,
                flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx,
                const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t lctx,
                const slong * perm, const ulong * shift, const ulong * stride);


/******************************************************************************

   Internal consistency checks

******************************************************************************/

void fmpz_mod_mpoly_remainder_strongtest(const fmpz_mod_mpoly_t r, const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
