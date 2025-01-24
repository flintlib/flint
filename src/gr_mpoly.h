/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_MPOLY_H
#define GR_MPOLY_H

#ifdef GR_MPOLY_INLINES_C
#define GR_MPOLY_INLINE
#else
#define GR_MPOLY_INLINE static inline
#endif

#include "mpoly_types.h"
#include "gr_vec.h"

#if FLINT_WANT_ASSERT
# include "mpoly.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    gr_ptr coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in ulong units */
    slong exps_alloc;       /* abs size in ulong units */
}
gr_mpoly_struct;

typedef gr_mpoly_struct gr_mpoly_t[1];

typedef struct
{
    gr_ctx_struct * cctx;
    mpoly_ctx_struct * mctx;
    char ** vars;
} _gr_mpoly_ctx_struct;

typedef gr_ctx_struct gr_mpoly_ctx_struct;
typedef gr_mpoly_ctx_struct gr_mpoly_ctx_t[1];

#define GR_MPOLY_MCTX(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->mctx)
#define GR_MPOLY_CCTX(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->cctx)
#define GR_MPOLY_VARS(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->vars)
#define GR_MPOLY_NVARS(ctx) (GR_MPOLY_MCTX(ctx)->nvars)

/* Context object */

void gr_mpoly_ctx_init_rand(gr_mpoly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring, slong max_nvars);
void gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx);

void gr_mpoly_ctx_init(gr_mpoly_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord);
void gr_mpoly_ctx_init_rand(gr_mpoly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring, slong max_nvars);

WARN_UNUSED_RESULT int gr_mpoly_ctx_set_gen_names(gr_mpoly_ctx_t ctx, const char ** s);
WARN_UNUSED_RESULT int gr_mpoly_gens(gr_vec_t res, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_gens_recursive(gr_vec_t vec, gr_mpoly_ctx_t ctx);

int gr_mpoly_ctx_write(gr_stream_t out, gr_mpoly_ctx_t ctx);
void gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx);

truth_t gr_mpoly_ctx_is_ring(gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_ctx_is_zero_ring(gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_ctx_is_commutative_ring(gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_ctx_is_integral_domain(gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_ctx_is_field(gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_ctx_is_threadsafe(gr_mpoly_ctx_t ctx);


/* Memory management */

GR_MPOLY_INLINE
void gr_mpoly_init(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx);
void gr_mpoly_init2( gr_mpoly_t A, slong alloc, gr_mpoly_ctx_t ctx);

GR_MPOLY_INLINE
void gr_mpoly_clear(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    _gr_vec_clear(A->coeffs, A->coeffs_alloc, GR_MPOLY_CCTX(ctx));

    if (A->coeffs_alloc > 0)
        flint_free(A->coeffs);
    if (A->exps_alloc > 0)
        flint_free(A->exps);
}

void _gr_mpoly_fit_length(
    gr_ptr * coeffs,
    slong * coeffs_alloc,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length,
    gr_ctx_t cctx);

void gr_mpoly_fit_length(gr_mpoly_t A, slong len, gr_mpoly_ctx_t ctx);
void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx);

void gr_mpoly_fit_length_fit_bits(
    gr_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    gr_mpoly_ctx_t ctx);

void gr_mpoly_fit_length_reset_bits(
    gr_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    gr_mpoly_ctx_t ctx);

/* todo: when to zero out coefficients? */
GR_MPOLY_INLINE
void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen, gr_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, GR_MPOLY_MCTX(ctx))*newlen <= A->exps_alloc);

    A->length = newlen;
}

GR_MPOLY_INLINE slong
gr_mpoly_length(const gr_mpoly_t x, gr_mpoly_ctx_t ctx)
{
    return x->length;
}

/* Basic manipulation */

GR_MPOLY_INLINE
void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, gr_mpoly_ctx_t ctx)
{
    FLINT_SWAP(gr_mpoly_struct, *A, *B);
}

GR_MPOLY_INLINE void
gr_mpoly_set_shallow(gr_mpoly_t res, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    *res = *poly;
}

WARN_UNUSED_RESULT int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx);

GR_MPOLY_INLINE WARN_UNUSED_RESULT
int gr_mpoly_zero(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    _gr_mpoly_set_length(A, 0, ctx);
    return GR_SUCCESS;
}

truth_t gr_mpoly_is_zero(const gr_mpoly_t A, gr_mpoly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mpoly_gen(gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx);
truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx);

truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx);

/* Container operations */

void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx);

void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx);

void gr_mpoly_sort_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_combine_like_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx);

truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx);
void gr_mpoly_assert_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx);

/* Random generation */

int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, gr_mpoly_ctx_t ctx);

GR_MPOLY_INLINE WARN_UNUSED_RESULT int
_gr_mpoly_randtest_default(gr_mpoly_t res, flint_rand_t state, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_randtest_bits(res, state, n_randint(state, 5), 1 + n_randint(state, 3), ctx);
}

/* Input and output */

int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, gr_mpoly_ctx_t ctx);
int gr_mpoly_write(gr_stream_t out, gr_mpoly_t poly, gr_mpoly_ctx_t ctx);
int gr_mpoly_print_pretty(const gr_mpoly_t A, gr_mpoly_ctx_t ctx);

/* Conversion */

WARN_UNUSED_RESULT int gr_mpoly_set_other(gr_mpoly_t res, gr_srcptr A, gr_ctx_t A_ctx, gr_mpoly_ctx_t ctx);

/* Constants */

WARN_UNUSED_RESULT int gr_mpoly_set_scalar(gr_mpoly_t A, gr_srcptr c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_ui(gr_mpoly_t A, ulong c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_si(gr_mpoly_t A, slong c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_fmpz(gr_mpoly_t A, const fmpz_t c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_fmpq(gr_mpoly_t A, const fmpq_t c, gr_mpoly_ctx_t ctx);

GR_MPOLY_INLINE WARN_UNUSED_RESULT
int gr_mpoly_one(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set_ui(A, 1, ctx);
}

truth_t gr_mpoly_is_one(const gr_mpoly_t A, gr_mpoly_ctx_t ctx);

/* Coefficient/exponent access */

WARN_UNUSED_RESULT int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, gr_mpoly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, gr_mpoly_ctx_t ctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mpoly_mul(gr_mpoly_t poly1, const gr_mpoly_t poly2, const gr_mpoly_t poly3, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_johnson(gr_mpoly_t poly1, const gr_mpoly_t poly2, const gr_mpoly_t poly3, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, gr_mpoly_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, gr_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
