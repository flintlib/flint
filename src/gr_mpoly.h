/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_MPOLY_H
#define GR_MPOLY_H

#ifdef GR_MPOLY_INLINES_C
#define GR_MPOLY_INLINE
#else
#define GR_MPOLY_INLINE static inline
#endif

#include "mpoly.h"
#include "gr_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    gr_ptr coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
}
gr_mpoly_struct;

typedef gr_mpoly_struct gr_mpoly_t[1];


/* Memory management */

GR_MPOLY_INLINE
void gr_mpoly_init(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits,
    const mpoly_ctx_t mctx, gr_ctx_t cctx);

void gr_mpoly_init2( gr_mpoly_t A, slong alloc,
    const mpoly_ctx_t mctx, gr_ctx_t cctx);

GR_MPOLY_INLINE
void gr_mpoly_clear(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    _gr_vec_clear(A->coeffs, A->coeffs_alloc, cctx);

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

void gr_mpoly_fit_length(gr_mpoly_t A, slong len, const mpoly_ctx_t mctx, gr_ctx_t cctx);

void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx);

void gr_mpoly_fit_length_fit_bits(
    gr_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx, gr_ctx_t cctx);

void gr_mpoly_fit_length_reset_bits(
    gr_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* todo: when to zero out coefficients? */
GR_MPOLY_INLINE
void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen,
                                               const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    FLINT_ASSERT(newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, mctx)*newlen <= A->exps_alloc);

    A->length = newlen;
}

/* Basic manipulation */

GR_MPOLY_INLINE
void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    FLINT_SWAP(gr_mpoly_struct, *A, *B);
}

WARN_UNUSED_RESULT int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx);

GR_MPOLY_INLINE WARN_UNUSED_RESULT
int gr_mpoly_zero(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    _gr_mpoly_set_length(A, 0, mctx, cctx);
    return GR_SUCCESS;
}

GR_MPOLY_INLINE
truth_t gr_mpoly_is_zero(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    if (A->length == 0)
        return T_TRUE;

    /* todo: skip when we have canonical representation */
    return _gr_vec_is_zero(A->coeffs, A->length, cctx);
}

WARN_UNUSED_RESULT int gr_mpoly_gen(gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx);
truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx);

truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* Container operations */

void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
void gr_mpoly_sort_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_combine_like_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx);
truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx);
void gr_mpoly_assert_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* Random generation */

int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* Input and output */

int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, const char ** x_in, const mpoly_ctx_t mctx, gr_ctx_t cctx);
int gr_mpoly_print_pretty(const gr_mpoly_t A, const char ** x_in, const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* Constants */

WARN_UNUSED_RESULT int gr_mpoly_set_scalar(gr_mpoly_t A, gr_srcptr c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_ui(gr_mpoly_t A, ulong c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_si(gr_mpoly_t A, slong c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_fmpz(gr_mpoly_t A, const fmpz_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_fmpq(gr_mpoly_t A, const fmpq_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx);

GR_MPOLY_INLINE WARN_UNUSED_RESULT
int gr_mpoly_one(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    return gr_mpoly_set_ui(A, 1, mctx, cctx);
}

/* todo: efficient version */
GR_MPOLY_INLINE
truth_t gr_mpoly_is_one(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    gr_mpoly_t t;
    truth_t res = T_UNKNOWN;
    gr_mpoly_init(t, mctx, cctx);
    if (gr_mpoly_one(t, mctx, cctx) == GR_SUCCESS)
        res = gr_mpoly_equal(A, t, mctx, cctx);
    gr_mpoly_clear(t, mctx, cctx);
    return res;
}

/* Coefficient/exponent access */

WARN_UNUSED_RESULT int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);

WARN_UNUSED_RESULT int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);

WARN_UNUSED_RESULT int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx);

/* Arithmetic */

WARN_UNUSED_RESULT int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx);

WARN_UNUSED_RESULT int gr_mpoly_mul(gr_mpoly_t poly1, const gr_mpoly_t poly2, const gr_mpoly_t poly3, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_johnson(gr_mpoly_t poly1, const gr_mpoly_t poly2, const gr_mpoly_t poly3, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx);

WARN_UNUSED_RESULT int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx);
WARN_UNUSED_RESULT int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx);

#ifdef __cplusplus
}
#endif

#endif
