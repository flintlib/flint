/*
    Copyright (C) 2017 - 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MPOLY_H
#define NMOD_MPOLY_H

#ifdef NMOD_MPOLY_INLINES_C
#define NMOD_MPOLY_INLINE
#else
#define NMOD_MPOLY_INLINE static inline
#endif

#include "thread_pool.h"
#include "mpoly.h"
#include "nmod_types.h"
#include "n_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

NMOD_MPOLY_INLINE
mp_limb_t * nmod_mpoly_term_coeff_ref(nmod_mpoly_t A, slong i,
                                                    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < A->length);
    return A->coeffs + i;
}


/* Internal type definitions *************************************************/

NMOD_MPOLY_INLINE n_poly_struct *
evil_cast_nmod_poly_to_n_poly(nmod_poly_struct * a)
{
    return (n_poly_struct *) a;
}

NMOD_MPOLY_INLINE const n_poly_struct *
evil_const_cast_nmod_poly_to_n_poly(const nmod_poly_struct * a)
{
    return (const n_poly_struct *) a;
}


/*
    nmod_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   nmod_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} nmod_mpoly_univar_struct;

typedef nmod_mpoly_univar_struct nmod_mpoly_univar_t[1];

/*
    nmod_mpolyu_t
    sparse univariates with nmod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   nmod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} nmod_mpolyu_struct;
typedef nmod_mpolyu_struct nmod_mpolyu_t[1];

/*
    nmod_mpolyn_t
    multivariates with n_poly_t coefficients
*/
typedef struct
{
   n_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} nmod_mpolyn_struct;
typedef nmod_mpolyn_struct nmod_mpolyn_t[1];

/*
    nmod_mpolyun_t
    sparse univariates with nmod_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    nmod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;   /* default bits to construct coeffs */
} nmod_mpolyun_struct;
typedef nmod_mpolyun_struct nmod_mpolyun_t[1];

/*
    nmod_mpolyd_t
    A dense mpoly is stored as a flat array of coefficients.
    Suppose deg_bounds = {r0, r1, r2}. The coefficient of the monomial with
    exponents {e0, e1, e2} is stored at the coefficient of index
        e2 + r2*(e1 + r1*(e0 + r0*0))
*/
typedef struct
{
    slong nvars;
    slong degb_alloc;
    slong * deg_bounds;
    slong coeff_alloc;
    mp_limb_t * coeffs;
} nmod_mpolyd_struct;

typedef nmod_mpolyd_struct nmod_mpolyd_t[1];

/* stack type used in gcd which are generally useful as well *****************/

typedef struct
{
    n_poly_struct ** poly_array;
    slong poly_alloc;
    slong poly_top;
    nmod_mpolyun_struct ** mpolyun_array;
    slong mpolyun_alloc;
    slong mpolyun_top;
    nmod_mpolyn_struct ** mpolyn_array;
    slong mpolyn_alloc;
    slong mpolyn_top;
    const nmod_mpoly_ctx_struct * ctx;
    flint_bitcnt_t bits;
} nmod_poly_stack_struct;

typedef nmod_poly_stack_struct nmod_poly_stack_t[1];

void nmod_poly_stack_init(nmod_poly_stack_t S, flint_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_poly_stack_clear(nmod_poly_stack_t S);

void nmod_poly_stack_set_ctx(nmod_poly_stack_t S,
                                                   const nmod_mpoly_ctx_t ctx);

n_poly_struct ** nmod_poly_stack_fit_request_poly(
                                                 nmod_poly_stack_t S, slong k);

nmod_mpolyun_struct ** nmod_poly_stack_fit_request_mpolyun(
                                                 nmod_poly_stack_t S, slong k);

nmod_mpolyn_struct ** nmod_poly_stack_fit_request_mpolyn(
                                                 nmod_poly_stack_t S, slong k);


NMOD_MPOLY_INLINE
n_poly_struct ** nmod_poly_stack_request_poly(nmod_poly_stack_t S, slong k)
{
    n_poly_struct ** poly_top;
    poly_top = nmod_poly_stack_fit_request_poly(S, k);
    S->poly_top += k;
    return poly_top;
}

NMOD_MPOLY_INLINE
n_poly_struct * nmod_poly_stack_take_top_poly(nmod_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    n_poly_struct ** poly_top;
    FLINT_ASSERT(S->poly_top + 1 <= S->poly_alloc);
    poly_top = S->poly_array + S->poly_top;
    S->poly_top += 1;
    return poly_top[0];
}

NMOD_MPOLY_INLINE
void nmod_poly_stack_give_back_poly(nmod_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->poly_top >= k);
    S->poly_top -= k;
}

NMOD_MPOLY_INLINE
slong nmod_poly_stack_size_poly(const nmod_poly_stack_t S)
{
    return S->poly_top;
}


NMOD_MPOLY_INLINE
nmod_mpolyun_struct ** nmod_poly_stack_request_mpolyun(nmod_poly_stack_t S, slong k)
{
    nmod_mpolyun_struct ** mpolyun_top;
    mpolyun_top = nmod_poly_stack_fit_request_mpolyun(S, k);
    S->mpolyun_top += k;
    return mpolyun_top;
}

NMOD_MPOLY_INLINE
nmod_mpolyun_struct * nmod_poly_stack_take_top_mpolyun(nmod_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    nmod_mpolyun_struct ** mpolyun_top;
    FLINT_ASSERT(S->mpolyun_top + 1 <= S->mpolyun_alloc);
    mpolyun_top = S->mpolyun_array + S->mpolyun_top;
    S->mpolyun_top += 1;
    return mpolyun_top[0];
}

NMOD_MPOLY_INLINE
void nmod_poly_stack_give_back_mpolyun(nmod_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->mpolyun_top >= k);
    S->mpolyun_top -= k;
}

NMOD_MPOLY_INLINE
slong nmod_poly_stack_size_mpolyun(const nmod_poly_stack_t S)
{
    return S->mpolyun_top;
}

NMOD_MPOLY_INLINE
nmod_mpolyn_struct ** nmod_poly_stack_request_mpolyn(nmod_poly_stack_t S, slong k)
{
    nmod_mpolyn_struct ** mpolyn_top;
    mpolyn_top = nmod_poly_stack_fit_request_mpolyn(S, k);
    S->mpolyn_top += k;
    return mpolyn_top;
}

NMOD_MPOLY_INLINE
nmod_mpolyn_struct * nmod_poly_stack_take_top_mpolyn(nmod_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    nmod_mpolyn_struct ** mpolyn_top;
    FLINT_ASSERT(S->mpolyn_top + 1 <= S->mpolyn_alloc);
    mpolyn_top = S->mpolyn_array + S->mpolyn_top;
    S->mpolyn_top += 1;
    return mpolyn_top[0];
}

NMOD_MPOLY_INLINE
void nmod_poly_stack_give_back_mpolyn(nmod_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->mpolyn_top >= k);
    S->mpolyn_top -= k;
}

NMOD_MPOLY_INLINE
slong nmod_poly_stack_size_mpolyn(const nmod_poly_stack_t S)
{
    return S->mpolyn_top;
}

/* Context object ************************************************************/

void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx,
                         slong nvars, const ordering_t ord, mp_limb_t modulus);

void nmod_mpoly_ctx_init_rand(nmod_mpoly_ctx_t ctx, flint_rand_t state,
                                           slong max_nvars, mp_limb_t modulus);

void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_ctx_set_modulus(nmod_mpoly_ctx_t ctx, mp_limb_t p)
{
    nmod_init(&ctx->mod, p);
}

NMOD_MPOLY_INLINE
slong nmod_mpoly_ctx_nvars(const nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

NMOD_MPOLY_INLINE
ordering_t nmod_mpoly_ctx_ord(const nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}

NMOD_MPOLY_INLINE
mp_limb_t nmod_mpoly_ctx_modulus(const nmod_mpoly_ctx_t ctx)
{
    return ctx->mod.n;
}

NMOD_MPOLY_INLINE
void nmod_mpoly_ctx_change_modulus(nmod_mpoly_ctx_t ctx, mp_limb_t modulus)
{
    nmod_init(&ctx->mod, modulus);
}

/*  Memory management ********************************************************/


NMOD_MPOLY_INLINE
void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

NMOD_MPOLY_INLINE
void nmod_mpoly_clear(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    if (A->coeffs_alloc > 0)
        flint_free(A->coeffs);

    if (A->exps_alloc > 0)
        flint_free(A->exps);
}

void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc,
                              flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_realloc(nmod_mpoly_t A,
                                      slong alloc, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_fit_length(nmod_mpoly_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_fit_length_fit_bits(nmod_mpoly_t A,
                   slong len, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_fit_length_reset_bits(nmod_mpoly_t A,
                   slong len, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void _nmod_mpoly_fit_length(
    mp_limb_t ** coeffs,
    slong * coeffs_alloc,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length)
{
    if (length > *coeffs_alloc)
    {
        *coeffs_alloc = FLINT_MAX(length, *coeffs_alloc*2);
        *coeffs = (mp_limb_t *) flint_realloc(*coeffs,
                                              *coeffs_alloc*sizeof(mp_limb_t));
    }

    if (N*length > *exps_alloc)
    {
        *exps_alloc = FLINT_MAX(N*length, *exps_alloc*2);
        *exps = (mp_limb_t *) flint_realloc(*exps, *exps_alloc*sizeof(ulong));
    }
}

NMOD_MPOLY_INLINE
void _nmod_mpoly_set_length(nmod_mpoly_t A, slong newlen,
                                                    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, ctx->minfo)*newlen <= A->exps_alloc);
    A->length = newlen;
}

NMOD_MPOLY_INLINE
void nmod_mpoly_truncate(nmod_mpoly_t A, slong newlen,
                                                   const nmod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        A->length = newlen;
    }
}


/* Input/output **************************************************************/

int nmod_mpoly_set_str_pretty(nmod_mpoly_t A, const char * str,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t A,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx);
#endif

int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx);

/*  Basic manipulation *******************************************************/

void nmod_mpoly_gen(nmod_mpoly_t A, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_is_gen(const nmod_mpoly_t A,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_equal(const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_swap(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    FLINT_SWAP(nmod_mpoly_struct, *A, *B);
}

/* Constants *****************************************************************/

int nmod_mpoly_is_ui(const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_ui(const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_ui(nmod_mpoly_t A, ulong c,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_fmpz(nmod_mpoly_t A, const fmpz_t c,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_zero(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
   _nmod_mpoly_set_length(A, 0, ctx);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_one(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_set_ui(A, UWORD(1), ctx);
}

int nmod_mpoly_equal_ui(const nmod_mpoly_t A,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
int nmod_mpoly_is_zero(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

NMOD_MPOLY_INLINE
int nmod_mpoly_is_one(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
   return nmod_mpoly_equal_ui(A, 1, ctx);
}

/* Degrees *******************************************************************/

NMOD_MPOLY_INLINE
int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_degrees_fmpz(fmpz ** degs, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t A, slong var,
                                                   const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

NMOD_MPOLY_INLINE
slong nmod_mpoly_degree_si(const nmod_mpoly_t A, slong var,
                                                   const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_total_degree_fits_si(const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_total_degree_fmpz(fmpz_t td, const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
slong nmod_mpoly_total_degree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_used_vars(int * used, const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}

/* Coefficients **************************************************************/

ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A,
                             const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_coeff_ui_monomial(nmod_mpoly_t A, ulong c,
                             const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_coeff_ui_ui(const nmod_mpoly_t A,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A,
                        ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A,
                      ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_coeff_ui_ui(nmod_mpoly_t A,
                       ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_coeff_vars_ui(nmod_mpoly_t C,
                 const nmod_mpoly_t A, const slong * vars, const ulong * exps,
                                     slong length, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE mp_limb_t nmod_mpoly_leadcoeff(
                                    nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs[0];
}

/* conversion ****************************************************************/

int nmod_mpoly_is_nmod_poly(const nmod_mpoly_t A,
                                        slong var, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_get_n_poly(n_poly_t A, const nmod_mpoly_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_get_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_set_nmod_poly(nmod_mpoly_t A, flint_bitcnt_t Abits,
                                        const mp_limb_t * Bcoeffs, slong Blen,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_n_poly_mod(nmod_mpoly_t A, const n_poly_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_nmod_poly(nmod_mpoly_t A, const nmod_poly_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

/* comparison ****************************************************************/

int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

int nmod_mpoly_is_canonical(const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
slong nmod_mpoly_length(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

void nmod_mpoly_resize(nmod_mpoly_t A, slong new_length,
                                                   const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_term_coeff_ui(const nmod_mpoly_t A, slong i,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
int nmod_mpoly_term_exp_fits_ui(const nmod_mpoly_t A, slong i,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_term_exp_fits_si(const nmod_mpoly_t A, slong i,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

void nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const nmod_mpoly_t A,
                                          slong i, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_term_exp_ui(ulong * exp, const nmod_mpoly_t A,
                                          slong i, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_term_exp_si(slong * exp, const nmod_mpoly_t A,
                                          slong i, const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_term_var_exp_ui(const nmod_mpoly_t A, slong i,
                                        slong var, const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_get_term_var_exp_si(const nmod_mpoly_t A, slong i,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_term_exp_fmpz(nmod_mpoly_t A, slong i,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_term_exp_ui(nmod_mpoly_t A, slong i,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_term(nmod_mpoly_t M, const nmod_mpoly_t A,
                                          slong i, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_term_monomial(nmod_mpoly_t M, const nmod_mpoly_t A,
                                          slong i, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_push_term_ui_ffmpz(nmod_mpoly_t A, ulong c,
                               const fmpz * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_push_term_ui_ui(nmod_mpoly_t A, ulong c,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_sort_terms(nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_combine_like_terms(nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_reverse(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_assert_canonical(const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_radix_sort1(nmod_mpoly_t A, slong left, slong right,
                              flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

void _nmod_mpoly_radix_sort(nmod_mpoly_t A, slong left, slong right,
                                    flint_bitcnt_t pos, slong N, ulong * cmpmask);

void _nmod_mpoly_push_exp_ffmpz(nmod_mpoly_t A,
                                 const fmpz * exp, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_push_exp_pfmpz(nmod_mpoly_t A,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_push_exp_ui(nmod_mpoly_t A,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state,
                 slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state,
                    slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state,
               slong length, flint_bitcnt_t exp_bits, const nmod_mpoly_ctx_t ctx);

ulong _nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly,
                                 const fmpz * exp, const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

ulong nmod_mpoly_get_term_ui_ui(const nmod_mpoly_t poly,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_max_degrees(ulong * max_degs, const ulong * exps,
                    slong len, slong bits, slong n, int deg, int rev, slong N);

void nmod_mpoly_max_degrees(ulong * max_degs,
                          const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_nmod(nmod_mpoly_t poly,
                                   const nmod_t c, const nmod_mpoly_ctx_t ctx);


ulong nmod_mpoly_get_coeff_ui(nmod_t x,
                 const nmod_mpoly_t poly, slong n, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_coeff_ui(nmod_mpoly_t poly,
                          slong n, ulong x, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_get_monomial(ulong * exps, const nmod_mpoly_t poly,
                                          slong n, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_set_monomial(nmod_mpoly_t poly,
                      slong n, const ulong * exps, const nmod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B,
                             const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_sub(nmod_mpoly_t A, const nmod_mpoly_t B,
                             const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

slong _nmod_mpoly_add(mp_limb_t * coeff1,       ulong * exp1,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                                  slong N, const ulong * cmpmask, nmod_t fctx);

slong _nmod_mpoly_sub(ulong * coeff1,       ulong * exp1,
                    const ulong * coeff2, const ulong * exp2, slong len2,
                    const ulong * coeff3, const ulong * exp3, slong len3,
                                  slong N, const ulong * cmpmask, nmod_t fctx);


/* Scalar operations *********************************************************/

void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A,
                const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_make_monic(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_scalar_mul_nmod_invertible(nmod_mpoly_t A,
                const nmod_mpoly_t B, mp_limb_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_scalar_addmul_ui(nmod_mpoly_t A,
                        const nmod_mpoly_t B, const nmod_mpoly_t C, ulong d,
                                                   const nmod_mpoly_ctx_t ctx);

/* Differention **************************************************************/

void nmod_mpoly_derivative(nmod_mpoly_t A,
                  const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

int _ff_poly_pow_fmpz_is_not_feasible(slong length, const fmpz_t e);

int _ff_poly_pow_ui_is_not_feasible(slong length, ulong e);

mp_limb_t _nmod_mpoly_eval_all_ui(const mp_limb_t * Acoeffs,
                 const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
                 const mp_limb_t * alphas, const mpoly_ctx_t mctx, nmod_t mod);

ulong nmod_mpoly_evaluate_all_ui(const nmod_mpoly_t A,
                               const ulong * vals, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                             slong var, ulong val, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_compose_nmod_poly(nmod_poly_t A,
                        const nmod_mpoly_t B, nmod_poly_struct * const * C,
                                                   const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_compose_mat(nmod_mpoly_t A,
                            const nmod_mpoly_t B, const fmpz_mat_t M,
                    const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC);

int nmod_mpoly_compose_nmod_mpoly_geobucket(nmod_mpoly_t A,
                    const nmod_mpoly_t B, nmod_mpoly_struct * const * C,
                    const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC);

int nmod_mpoly_compose_nmod_mpoly_horner(nmod_mpoly_t A,
                    const nmod_mpoly_t B, nmod_mpoly_struct * const * C,
                    const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC);

int nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A,
                    const nmod_mpoly_t B, nmod_mpoly_struct * const * C,
                    const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC);

void nmod_mpoly_compose_nmod_mpoly_gen(nmod_mpoly_t A,
                    const nmod_mpoly_t B, const slong * c,
                    const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC);


/* Multiplication ************************************************************/

void nmod_mpoly_mul(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_mul_johnson(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_mul_array(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_mul_array_threaded(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_mul_dense(nmod_mpoly_t A,
       const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx);

slong _nmod_mpoly_mul_johnson(nmod_mpoly_t A,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
             flint_bitcnt_t bits, slong N, const ulong * cmpmask, nmod_t fctx);

void _nmod_mpoly_mul_johnson_maxfields(nmod_mpoly_t A,
                                 const nmod_mpoly_t B, fmpz * maxBfields,
                                 const nmod_mpoly_t C, fmpz * maxCfields,
                                                   const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_mul_heap_threaded_pool_maxfields(nmod_mpoly_t A,
           const nmod_mpoly_t B, fmpz * maxBfields,
           const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

int _nmod_mpoly_mul_array_DEG(nmod_mpoly_t A,
                                 const nmod_mpoly_t B, fmpz * maxBfields,
                                 const nmod_mpoly_t C, fmpz * maxCfields,
                                                   const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_mul_array_LEX(nmod_mpoly_t A,
                                 const nmod_mpoly_t B, fmpz * maxBfields,
                                 const nmod_mpoly_t C, fmpz * maxCfields,
                                                   const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_mul_array_threaded_pool_DEG(nmod_mpoly_t A,
           const nmod_mpoly_t B, fmpz * maxBfields,
           const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

int _nmod_mpoly_mul_array_threaded_pool_LEX(nmod_mpoly_t A,
           const nmod_mpoly_t B, fmpz * maxBfields,
           const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

int _nmod_mpoly_mul_dense(nmod_mpoly_t P,
                                 const nmod_mpoly_t A, fmpz * maxAfields,
                                 const nmod_mpoly_t B, fmpz * maxBfields,
                                                   const nmod_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

void _nmod_mpoly_pow_rmul(nmod_mpoly_t A, const mp_limb_t * Bcoeffs,
                            const ulong * Bexps, slong Blen, ulong k, slong N,
                            const ulong * cmpmask, nmod_t mod, nmod_mpoly_t T);

void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B,
                                          ulong k, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B,
                                   const fmpz_t k, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                          ulong k, const nmod_mpoly_ctx_t ctx);


/* Division ******************************************************************/

int nmod_mpoly_divides(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_divides_threaded_pool(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

#if FLINT_KNOW_STRONG_ORDER
#define nmod_mpoly_divides_heap_threaded nmod_mpoly_divides_heap_threaded
int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

#define _nmod_mpoly_divides_heap_threaded_pool _nmod_mpoly_divides_heap_threaded_pool
int _nmod_mpoly_divides_heap_threaded_pool(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);
#endif

int nmod_mpoly_divides_dense(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_div(nmod_mpoly_t Q,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_divrem(nmod_mpoly_t Q, nmod_mpoly_t R,
       const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R,
                const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_divexact(nmod_mpoly_t Q, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_divides(Q, A, B, ctx))
        return;

    flint_throw(FLINT_ERROR, "nmod_mpoly_divexact: nonexact division");
}

int _nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q,
                    const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                    const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
                                                                  nmod_t fctx);

void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t Q,
                                 const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

void
nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** Q, nmod_mpoly_t R,
    const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len,
                                                   const nmod_mpoly_ctx_t ctx);

/* Square root ***************************************************************/

int nmod_mpoly_sqrt_heap(nmod_mpoly_t Q, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
int nmod_mpoly_sqrt(nmod_mpoly_t Q, const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_sqrt_heap(Q, A, ctx);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_is_square(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    int res;
    nmod_mpoly_t Q;
    nmod_mpoly_init(Q, ctx);
    res = nmod_mpoly_sqrt_heap(Q, A, ctx);
    nmod_mpoly_clear(Q, ctx);
    return res;
}

int nmod_mpoly_quadratic_root(nmod_mpoly_t Q, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

void nmod_mpoly_term_content(nmod_mpoly_t M, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_content_vars(nmod_mpoly_t g, const nmod_mpoly_t A,
                  slong * vars, slong vars_length, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int _nmod_mpoly_gcd_algo_small(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar, /* could be NULL */
    nmod_mpoly_t Bbar, /* could be NULL */
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo);

int _nmod_mpoly_gcd_algo(nmod_mpoly_t G, nmod_mpoly_t Abar,
                nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B,
                                const nmod_mpoly_ctx_t ctx, unsigned int algo);

int nmod_mpoly_gcd_cofactors(nmod_mpoly_t G,
                nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd_hensel(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd_zippel2(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_deflation(fmpz * shift, fmpz * stride,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_deflate(nmod_mpoly_t A, const nmod_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_inflate(nmod_mpoly_t A, const nmod_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const nmod_mpoly_ctx_t ctx);


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

void mpoly_void_ring_init_nmod_mpoly_ctx(mpoly_void_ring_t R,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyl_lead_coeff(nmod_mpoly_t c, const nmod_mpoly_t A,
                                   slong num_vars, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyl_content(nmod_mpoly_t g, const nmod_mpoly_t A,
                                   slong num_vars, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_to_nmod_poly_deflate(nmod_poly_t A, const nmod_mpoly_t B,
                        slong var, const ulong * Bshift, const ulong * Bstride,
                                                   const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_from_nmod_poly_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits,
                         const nmod_poly_t B, slong var, const ulong * Ashift,
                            const ulong * Astride, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_repack_bits(nmod_mpoly_t A, const nmod_mpoly_t B,
                             flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_repack_bits_inplace(nmod_mpoly_t A,
                             flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx);

typedef struct
{
    slong nvars;
    slong * perm;
} nmod_mpolyd_ctx_struct;

typedef nmod_mpolyd_ctx_struct nmod_mpolyd_ctx_t[1];


/* data is passed to the threaded mul/div functions via a stripe struct */

typedef struct _nmod_mpoly_stripe_struct
{
    char * big_mem;
    slong big_mem_alloc;
    const nmod_mpoly_ctx_struct * ctx;
    slong N;
    flint_bitcnt_t bits;
    nmod_t mod;
    mp_limb_t lc_minus_inv;
    const ulong * cmpmask;
    slong * startidx;
    slong * endidx;
    ulong * emin;
    ulong * emax;
    int upperclosed;
} nmod_mpoly_stripe_struct;

typedef nmod_mpoly_stripe_struct nmod_mpoly_stripe_t[1];


/* Univariates ***************************************************************/

void nmod_mpoly_univar_init(nmod_mpoly_univar_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_univar_clear(nmod_mpoly_univar_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_univar_fit_length(nmod_mpoly_univar_t A,
                                     slong length, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_univar_print_pretty(const nmod_mpoly_univar_t A,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_univar_assert_canonical(nmod_mpoly_univar_t A,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_univar_zero(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

void nmod_mpoly_univar_set_coeff_ui(nmod_mpoly_univar_t A,
                    ulong e, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_to_univar(nmod_mpoly_univar_t A,
                  const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_from_univar(nmod_mpoly_t A, flint_bitcnt_t Abits,
           const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_from_univar(nmod_mpoly_t A,
           const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_univar_swap(nmod_mpoly_univar_t A, nmod_mpoly_univar_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_univar_struct t = *A;
    *A = *B;
    *B = t;
}

NMOD_MPOLY_INLINE
int nmod_mpoly_univar_degree_fits_si(const nmod_mpoly_univar_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

NMOD_MPOLY_INLINE
slong nmod_mpoly_univar_length(const nmod_mpoly_univar_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

NMOD_MPOLY_INLINE
slong nmod_mpoly_univar_get_term_exp_si(nmod_mpoly_univar_t A, slong i,
                                                    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_univar_get_term_coeff(nmod_mpoly_t c,
              const nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    nmod_mpoly_set(c, A->coeffs + i, ctx);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_univar_swap_term_coeff(nmod_mpoly_t c,
                    nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    nmod_mpoly_swap(c, A->coeffs + i, ctx);
}

int nmod_mpoly_univar_pseudo_gcd(nmod_mpoly_univar_t Gx,
        const nmod_mpoly_univar_t Ax, const nmod_mpoly_univar_t Bx,
                                               const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_univar_resultant(nmod_mpoly_t R,
        const nmod_mpoly_univar_t Ax, const nmod_mpoly_univar_t Bx,
                                               const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_univar_discriminant(nmod_mpoly_t D,
             const nmod_mpoly_univar_t Fx, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_resultant(nmod_mpoly_t R,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                    slong var, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_discriminant(nmod_mpoly_t R,
          const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx);

/* Helpers for array methods *************************************************/

void _nmod_mpoly_mul_array_chunked_LEX(nmod_mpoly_t P,
                             const nmod_mpoly_t A, const nmod_mpoly_t B,
                              const ulong * mults, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_mul_array_chunked_DEG(nmod_mpoly_t P,
                             const nmod_mpoly_t A, const nmod_mpoly_t B,
                                       ulong degb, const nmod_mpoly_ctx_t ctx);

void _nmod_mpoly_addmul_array1_ulong1(ulong * poly1,
                          const ulong * poly2, const ulong * exp2, slong len2,
                          const ulong * poly3, const ulong * exp3, slong len3);

void _nmod_mpoly_addmul_array1_ulong2(ulong * poly1,
                          const ulong * poly2, const ulong * exp2, slong len2,
                          const ulong * poly3, const ulong * exp3, slong len3);

void _nmod_mpoly_addmul_array1_ulong3(ulong * poly1,
                          const ulong * poly2, const ulong * exp2, slong len2,
                          const ulong * poly3, const ulong * exp3, slong len3);

slong nmod_mpoly_append_array_sm1_LEX(nmod_mpoly_t P, slong Plen,
         ulong * coeff_array, const ulong * mults, slong num, slong array_size,
                                        slong top, const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm2_LEX(nmod_mpoly_t P, slong Plen,
         ulong * coeff_array, const ulong * mults, slong num, slong array_size,
                                        slong top, const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm3_LEX(nmod_mpoly_t P, slong Plen,
         ulong * coeff_array, const ulong * mults, slong num, slong array_size,
                                        slong top, const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm1_DEGLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm2_DEGLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm3_DEGLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm1_DEGREVLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm2_DEGREVLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

slong nmod_mpoly_append_array_sm3_DEGREVLEX(nmod_mpoly_t P, slong Plen,
                     ulong * coeff_array, slong top, slong nvars, slong degb,
                                                   const nmod_mpoly_ctx_t ctx);

/* mpolyd ********************************************************************/

void nmod_mpolyd_ctx_init(nmod_mpolyd_ctx_t dctx, slong nvars);

void nmod_mpolyd_ctx_clear(nmod_mpolyd_ctx_t dctx);

NMOD_MPOLY_INLINE void nmod_mpolyd_swap(nmod_mpolyd_t poly1,
                                nmod_mpolyd_t poly2)
{
   nmod_mpolyd_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

int nmod_mpolyd_set_degbounds(nmod_mpolyd_t A, slong * bounds);

int nmod_mpolyd_set_degbounds_perm(nmod_mpolyd_t A,
                                 const nmod_mpolyd_ctx_t dctx, slong * bounds);

void nmod_mpoly_convert_to_nmod_mpolyd(
                                  nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_convert_to_nmod_mpolyd_degbound(
                                  nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_convert_from_nmod_mpolyd(
                                 nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx,
                          const nmod_mpolyd_t B, const nmod_mpolyd_ctx_t dctx);

void nmod_mpolyd_init(nmod_mpolyd_t poly, slong nvars);

void nmod_mpolyd_fit_length(nmod_mpolyd_t poly, slong len);

void nmod_mpolyd_zero(nmod_mpolyd_t poly);

void nmod_mpolyd_set_nvars(nmod_mpolyd_t poly, slong nvars);

void nmod_mpolyd_set(nmod_mpolyd_t A, const nmod_mpolyd_t B);

void nmod_mpolyd_clear(nmod_mpolyd_t poly);

void nmod_mpolyd_print(nmod_mpolyd_t poly);

slong nmod_mpolyd_length(const nmod_mpolyd_t A);


/* mpolyu ********************************************************************/

void nmod_mpolyu_init(nmod_mpolyu_t A, flint_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                  const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length,
                                                  const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_degrees_si(
    slong * degs,
    const nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_repack_bits_inplace(
    nmod_mpolyu_t A,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx);

nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A,
                                       ulong pow, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_shift_right(nmod_mpolyu_t A, ulong s);

void nmod_mpolyu_shift_left(nmod_mpolyu_t A, ulong s);

int nmod_mpolyu_content_mpoly(nmod_mpoly_t g,
                           const nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c,
                                                  const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B,
                                                  const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_cvtfrom_poly(nmod_mpolyu_t A, nmod_poly_t a,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(
                nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx,
                const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride,
                       const thread_pool_handle * handles, slong num_handles);

void nmod_mpoly_from_mpolyu_perm_inflate(
            nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx,
                        const nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx,
                const slong * perm, const ulong * shift, const ulong * stride);

int nmod_mpolyuu_divides(nmod_mpolyu_t Q, const nmod_mpolyu_t A,
           const nmod_mpolyu_t B, slong nmainvars, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_divexact_mpoly_inplace(nmod_mpolyu_t A,
                                   nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_mul_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                   nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_mul_mpoly_inplace(nmod_mpolyu_t A, nmod_mpoly_t c,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyu_gcdm_zippel(nmod_mpolyu_t G, nmod_mpolyu_t Abar,
                       nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B,
                                 nmod_mpoly_ctx_t ctx, flint_rand_t randstate);

NMOD_MPOLY_INLINE mp_limb_t nmod_mpolyu_leadcoeff(
                                   nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpoly_leadcoeff(A->coeffs + 0, ctx);
}

/* mpolyn ********************************************************************/

void nmod_mpolyn_init(nmod_mpolyn_t A, flint_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_clear(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_swap(nmod_mpolyn_t A, nmod_mpolyn_t B);

void nmod_mpolyn_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_is_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_print_pretty(const nmod_mpolyn_t A, const char ** x_in,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_fit_length(nmod_mpolyn_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_set_length(nmod_mpolyn_t A, slong newlen,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_fit_bits(nmod_mpolyn_t A, slong bits,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_is_canonical(const nmod_mpolyn_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_set(nmod_mpolyn_t A, const nmod_mpolyn_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_set_mpoly(nmod_mpolyn_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_cvtfrom_mpolyn(nmod_mpoly_t A, const nmod_mpolyn_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_mul_poly(nmod_mpolyn_t A, const nmod_mpolyn_t B,
                                 const n_poly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_cvtto_mpolyn(nmod_mpolyn_t A, const nmod_mpoly_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE mp_limb_t nmod_mpolyn_leadcoeff(nmod_mpolyn_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    n_poly_struct * leadpoly;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_poly_degree(A->coeffs + 0) >= 0);

    leadpoly = A->coeffs + 0;
    return leadpoly->coeffs[leadpoly->length - 1];
}

NMOD_MPOLY_INLINE n_poly_struct * nmod_mpolyn_leadcoeff_poly(
                                   nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* mpolyun *******************************************************************/

void nmod_mpolyun_init(nmod_mpolyun_t A, flint_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B);

void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyun_is_canonical(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s);

void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s);

slong nmod_mpolyn_lastdeg(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

slong nmod_mpolyun_lastdeg(nmod_mpolyun_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_one(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_one(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

mp_limb_t nmod_mpolyun_leadcoeff_last(nmod_mpolyun_t A,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_set_mod(nmod_mpolyn_t A, const nmod_t mod);

void nmod_mpolyun_set_mod(nmod_mpolyun_t A, const nmod_t mod);

int nmod_mpolyn_is_nonzero_nmod(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyun_is_nonzero_nmod(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_scalar_mul_nmod(
    nmod_mpolyn_t A,
    mp_limb_t c,
    const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, mp_limb_t c,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_mul_last(nmod_mpolyn_t A, n_poly_t b,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_mul_last(nmod_mpolyun_t A, n_poly_t b,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_equal(
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx);

int nmod_mpolyun_equal(const nmod_mpolyun_t A,
                           const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_cvtto_mpolyun(nmod_mpolyun_t A, const nmod_mpolyu_t B,
                                          slong k, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_cvtfrom_mpolyun(nmod_mpolyu_t A, const nmod_mpolyun_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_mul_poly(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                                 const n_poly_t c, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_content_last(n_poly_t a, nmod_mpolyn_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_content_last(n_poly_t a, nmod_mpolyun_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_divexact_last(nmod_mpolyn_t A, n_poly_t b,
                                                    const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, n_poly_t b,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_divides(nmod_mpolyn_t Q, const nmod_mpolyn_t A,
                            const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx);

#if FLINT_KNOW_STRONG_ORDER
#define nmod_mpolyn_divides_threaded_pool nmod_mpolyn_divides_threaded_pool
int nmod_mpolyn_divides_threaded_pool(nmod_mpolyn_t Q,
    const nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);
#endif

int nmod_mpolyun_divides(nmod_mpolyun_t Q, const nmod_mpolyun_t A,
                           const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_to_mpolyun_perm_deflate_threaded_pool(
                nmod_mpolyun_t A, const nmod_mpoly_ctx_t uctx,
                const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride,
                        const thread_pool_handle * handles, slong num_handles);

void nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(nmod_mpolyn_t A,
  const nmod_mpoly_ctx_t nctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                 const slong * perm, const ulong * shift, const ulong * stride,
                        const thread_pool_handle * handles, slong num_handles);

void nmod_mpoly_from_mpolyun_perm_inflate(nmod_mpoly_t A,
     flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx, const nmod_mpolyun_t B,
        const nmod_mpoly_ctx_t uctx, const slong * perm, const ulong * shift,
                                                        const ulong * stride);

void nmod_mpoly_from_mpolyn_perm_inflate(nmod_mpoly_t A,
                        flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx,
                        const nmod_mpolyn_t B, const nmod_mpoly_ctx_t nctx,
                const slong * perm, const ulong * shift, const ulong * stride);

NMOD_MPOLY_INLINE mp_limb_t nmod_mpolyun_leadcoeff(
                                  nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff(A->coeffs + 0, ctx);
}

NMOD_MPOLY_INLINE n_poly_struct * nmod_mpolyun_leadcoeff_poly(
                                  nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff_poly(A->coeffs + 0, ctx);
}


/* GCD ***********************************************************************/

int mpoly_gcd_get_use_first(slong rGdeg, slong Adeg, slong Bdeg,
                                                               slong gammadeg);

int nmod_mpoly_gcd_get_use_new(slong rGdeg, slong Adeg, slong Bdeg,
             slong gammadeg, slong degxAB, slong degyAB, slong numABgamma,
             const n_polyun_t G, const n_polyun_t Abar, const n_polyun_t Bbar);

void nmod_mpolyu_setform_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B,
                                                   const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_gcd_brown_smprime_bivar(
               nmod_mpolyn_t G, nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar,
               nmod_mpolyn_t A, nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx,
                                                         nmod_poly_stack_t Sp);

int nmod_mpolyn_gcd_brown_smprime(nmod_mpolyn_t G,
                                  nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar,
                                 nmod_mpolyn_t A, nmod_mpolyn_t B, slong var,
                         const nmod_mpoly_ctx_t ctx, const mpoly_gcd_info_t I,
                                                         nmod_poly_stack_t Sp);

int nmod_mpolyn_gcd_brown_smprime_threaded_pool(nmod_mpolyn_t G,
                                nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar,
                               nmod_mpolyn_t A, nmod_mpolyn_t B, slong var,
                         const nmod_mpoly_ctx_t ctx, const mpoly_gcd_info_t I,
                        const thread_pool_handle * handles, slong num_workers);

int nmod_mpolyn_gcd_brown_lgprime(nmod_mpolyn_t G,
                                 nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar,
                                 nmod_mpolyn_t A, nmod_mpolyn_t B, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

typedef enum {
    nmod_gcds_success,
    nmod_gcds_form_main_degree_too_high,
    nmod_gcds_form_wrong,
    nmod_gcds_no_solution,
    nmod_gcds_scales_not_found,
    nmod_gcds_eval_point_not_found,
    nmod_gcds_eval_gcd_deg_too_high
} nmod_gcds_ret_t;

nmod_gcds_ret_t nmod_mpolyu_gcds_zippel(nmod_mpolyu_t G,
                           nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpolyu_t f,
                                        slong var, const nmod_mpoly_ctx_t ctx,
                                     flint_rand_t randstate, slong * degbound);

int nmod_mpolyu_gcdp_zippel(nmod_mpolyu_t G, nmod_mpolyu_t Abar,
             nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B, slong var,
                           const nmod_mpoly_ctx_t ctx, flint_rand_t randstate);

void nmod_mpoly_to_mpolyl_perm_deflate(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t lctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

void nmod_mpoly_from_mpolyl_perm_inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

int nmod_mpolyl_gcd_zippel_smprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t ctx);

int nmod_mpolyl_gcd_zippel_lgprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t ctx);

int nmod_mpolyl_gcd_hensel_smprime(
    nmod_mpoly_t G, slong Gdeg,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx);

int nmod_mpolyl_gcd_hensel_medprime(
    nmod_mpoly_t G, slong Gdeg,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t smA,
    const nmod_mpoly_t smB,
    const nmod_mpoly_ctx_t smctx);

void _nmod_mpoly_monomial_evals_cache(n_poly_t E,
                    const ulong * Aexps, flint_bitcnt_t Abits, slong Alen,
                    const mp_limb_t * betas, slong start, slong stop,
                                          const mpoly_ctx_t mctx, nmod_t mod);

void _nmod_mpoly_monomial_evals2_cache(n_polyun_t E,
          const ulong * Aexps, flint_bitcnt_t Abits, slong Alen,
          const mp_limb_t * betas, slong m, const mpoly_ctx_t ctx, nmod_t mod);

/* interp ********************************************************************/

void _nmod_poly_eval2_pow(mp_limb_t * vp, mp_limb_t * vm,
                                   n_poly_t P, n_poly_t alphapow, nmod_t fctx);

void nmod_mpolyn_interp_reduce_2sm_poly(n_poly_t E,
                         n_poly_t F, const nmod_mpolyn_t A, n_poly_t alphapow,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_2sm_poly(slong * lastdeg_,
                       nmod_mpolyn_t F, const n_poly_t A, const n_poly_t B,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_interp_crt_2sm_poly(slong * lastdeg_,
                     nmod_mpolyn_t F, nmod_mpolyn_t T, const n_poly_t A,
                             const n_poly_t B, const n_poly_t modulus,
                               n_poly_t alphapow,  const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_sm_bpoly(nmod_mpolyn_t F,
                                      n_bpoly_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_interp_crt_sm_bpoly(slong * lastdeg,
            nmod_mpolyn_t F, nmod_mpolyn_t T, n_bpoly_t A, n_poly_t modulus,
                                n_poly_t alphapow, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_reduce_2sm_mpolyn(nmod_mpolyn_t E,
                            nmod_mpolyn_t F, nmod_mpolyn_t A, slong var,
                                n_poly_t alphapow, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_2sm_mpolyn(slong * lastdeg,
                      nmod_mpolyn_t T, nmod_mpolyn_t A, nmod_mpolyn_t B,
                       slong var, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_interp_crt_2sm_mpolyn(slong * lastdeg,
                        nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_mpolyn_t A,
                        nmod_mpolyn_t B, slong var, n_poly_t modulus,
                               n_poly_t alphapow, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_interp_reduce_sm_mpolyu(nmod_mpolyu_t B,
               nmod_mpolyun_t A, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_sm_mpoly(nmod_mpolyn_t A,
                            const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_interp_lift_sm_mpolyu(nmod_mpolyun_t A,
                           const nmod_mpolyu_t B, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_interp_crt_sm_mpoly(slong * lastdeg, nmod_mpolyn_t F,
                            nmod_mpolyn_t T, nmod_mpoly_t A, n_poly_t modulus,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyun_interp_crt_sm_mpolyu(slong * lastdeg,
             nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_mpolyu_t A,
                n_poly_t modulus, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyn_interp_mcrt_sm_mpoly(slong * lastdeg_,
                nmod_mpolyn_t F, const nmod_mpoly_t A, const n_poly_t modulus,
                                n_poly_t alphapow, const nmod_mpoly_ctx_t ctx);

/* geobuckets ****************************************************************/

typedef struct nmod_mpoly_geobucket
{
    nmod_mpoly_struct polys[FLINT_BITS/2];
    nmod_mpoly_struct temps[FLINT_BITS/2];
    slong length;
} nmod_mpoly_geobucket_struct;

typedef nmod_mpoly_geobucket_struct nmod_mpoly_geobucket_t[1];

void nmod_mpoly_geobucket_init(nmod_mpoly_geobucket_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_clear(nmod_mpoly_geobucket_t B,
                                                   const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_empty(nmod_mpoly_t p,
                         nmod_mpoly_geobucket_t B, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_fit_length(nmod_mpoly_geobucket_t B,
                                          slong i, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_set(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_add(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_geobucket_sub(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal consistency checks

******************************************************************************/

void nmod_mpoly_remainder_strongtest(const nmod_mpoly_t r, const nmod_mpoly_t g, const nmod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

