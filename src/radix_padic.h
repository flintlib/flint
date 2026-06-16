/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef RADIX_PADIC_H
#define RADIX_PADIC_H

#ifdef RADIX_PADIC_INLINES_C
#define RADIX_PADIC_INLINE
#else
#define RADIX_PADIC_INLINE static inline
#endif

#include "radix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RADIX_PADIC_EXACT     WORD_MAX            /* N value marking an exact element */
#define RADIX_PADIC_PREC_INF  WORD_MAX            /* context precision: "infinite" */
#define RADIX_PADIC_ERR_MAX   (WORD_MAX / 4)      /* largest finite precision we store */

/* context option flags */
#define RADIX_PADIC_SIGNED    1                   /* allow signed units for exact elements */
#define RADIX_PADIC_NO_ERROR  2                   /* floating-point mode: never track error (future) */
#define RADIX_PADIC_DECIMAL   4                   /* print the unit as a decimal integer */
#define RADIX_PADIC_TEST_LIMITS 8                 /* randtest: exercise the v/N representability limits */

typedef struct
{
    radix_integer_struct u;     /* unit (limbs in radix B = p^e) */
    slong v;                    /* valuation, in powers of p */
    slong N;                    /* O(p^N); RADIX_PADIC_EXACT if exact */
}
radix_padic_struct;

typedef radix_padic_struct radix_padic_t[1];

#define RADIX_PADIC_UNIT(x)  (&((x)->u))
#define RADIX_PADIC_VAL(x)   ((x)->v)
#define RADIX_PADIC_N(x)     ((x)->N)

typedef struct
{
    radix_struct radix;     /* digit radix p, limb radix p^e */
    ulong p;                /* the prime */
    slong prec_abs;         /* default absolute precision (RADIX_PADIC_PREC_INF for none) */
    slong prec_rel;         /* default relative precision (RADIX_PADIC_PREC_INF for none) */
    int flags;
}
radix_padic_ctx_struct;

typedef radix_padic_ctx_struct radix_padic_ctx_t[1];

#define GR_RADIX_PADIC_CTX(ctx)  ((radix_padic_ctx_struct *) (GR_CTX_DATA_AS_PTR(ctx)))
#define RADIX_PADIC_CTX_RADIX(ctx)  (&(GR_RADIX_PADIC_CTX(ctx)->radix))
#define RADIX_PADIC_CTX_PREC_ABS(ctx) (GR_RADIX_PADIC_CTX(ctx)->prec_abs)
#define RADIX_PADIC_CTX_PREC_REL(ctx) (GR_RADIX_PADIC_CTX(ctx)->prec_rel)
#define RADIX_PADIC_CTX_FLAGS(ctx)  (GR_RADIX_PADIC_CTX(ctx)->flags)
#define RADIX_PADIC_CTX_SIGNED(ctx) (GR_RADIX_PADIC_CTX(ctx)->flags & RADIX_PADIC_SIGNED)
#define RADIX_PADIC_CTX_DECIMAL(ctx) (GR_RADIX_PADIC_CTX(ctx)->flags & RADIX_PADIC_DECIMAL)

/* Context management */
int gr_ctx_init_radix_padic(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags);
void radix_padic_ctx_clear(gr_ctx_t ctx);
int radix_padic_ctx_write(gr_stream_t out, gr_ctx_t ctx);

/* Memory management */
void radix_padic_init(radix_padic_t res, gr_ctx_t ctx);
void radix_padic_clear(radix_padic_t res, gr_ctx_t ctx);
void radix_padic_swap(radix_padic_t x, radix_padic_t y, gr_ctx_t ctx);
void radix_padic_set_shallow(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);

/* Error term access */
slong radix_padic_get_error(const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_is_exact(const radix_padic_t x, gr_ctx_t ctx);

int _radix_padic_finalize(radix_padic_t res, gr_ctx_t ctx);

/* Basic assignment and predicates */
int radix_padic_randtest(radix_padic_t res, flint_rand_t state, gr_ctx_t ctx);
int radix_padic_write(gr_stream_t out, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_zero(radix_padic_t res, gr_ctx_t ctx);
int radix_padic_one(radix_padic_t res, gr_ctx_t ctx);
int radix_padic_set(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_set_ui(radix_padic_t res, ulong x, gr_ctx_t ctx);
int radix_padic_set_si(radix_padic_t res, slong x, gr_ctx_t ctx);
int radix_padic_set_fmpz(radix_padic_t res, const fmpz_t x, gr_ctx_t ctx);
int radix_padic_exact_set_ui(radix_padic_t res, ulong x, gr_ctx_t ctx);
int radix_padic_exact_set_si(radix_padic_t res, slong x, gr_ctx_t ctx);
int radix_padic_exact_set_fmpz(radix_padic_t res, const fmpz_t x, gr_ctx_t ctx);
int radix_padic_get_fmpz(fmpz_t res, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_neg(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_add(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);
int radix_padic_sub(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);
int radix_padic_mul(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);
int radix_padic_inv(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_div(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);
int radix_padic_sqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);
int radix_padic_rsqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_is_square(const radix_padic_t x, gr_ctx_t ctx);

truth_t radix_padic_is_zero(const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_is_one(const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_is_neg_one(const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_is_invertible(const radix_padic_t x, gr_ctx_t ctx);
truth_t radix_padic_equal(const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);

int radix_padic_dot(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1,
    const radix_padic_struct * vec2, slong len, gr_ctx_t ctx);
int radix_padic_dot_rev(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1,
    const radix_padic_struct * vec2, slong len, gr_ctx_t ctx);

int radix_padic_dot_strided(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);
int radix_padic_dot_strided_delayed(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);
int radix_padic_dot_strided_naive(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);

/* For testing */
int _radix_padic_add_sub_reference(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, int sub, gr_ctx_t ctx);
int _radix_padic_mul_reference(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
