/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_RADIX_H
#define PADIC_RADIX_H

#ifdef PADIC_RADIX_INLINES_C
#define PADIC_RADIX_INLINE
#else
#define PADIC_RADIX_INLINE static inline
#endif

#include "radix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PADIC_RADIX_EXACT     WORD_MAX            /* N value marking an exact element */
#define PADIC_RADIX_PREC_INF  WORD_MAX            /* context precision: "infinite" */
#define PADIC_RADIX_ERR_MAX   (WORD_MAX / 4)      /* largest finite precision we store */

/* context option flags */
#define PADIC_RADIX_SIGNED    1                   /* allow signed units for exact elements */
#define PADIC_RADIX_NO_ERROR  2                   /* floating-point mode: never track error (future) */
#define PADIC_RADIX_DECIMAL   4                   /* print the unit as a decimal integer */
#define PADIC_RADIX_TEST_LIMITS 8                 /* randtest: exercise the v/N representability limits */

typedef struct
{
    radix_integer_struct u;     /* unit (limbs in radix B = p^e) */
    slong v;                    /* valuation, in powers of p */
    slong N;                    /* O(p^N); PADIC_RADIX_EXACT if exact */
}
padic_radix_struct;

typedef padic_radix_struct padic_radix_t[1];

#define PADIC_RADIX_UNIT(x)  (&((x)->u))
#define PADIC_RADIX_VAL(x)   ((x)->v)
#define PADIC_RADIX_N(x)     ((x)->N)

typedef struct
{
    radix_struct radix;     /* digit radix p, limb radix p^e */
    ulong p;                /* the prime */
    slong prec_abs;         /* default absolute precision (PADIC_RADIX_PREC_INF for none) */
    slong prec_rel;         /* default relative precision (PADIC_RADIX_PREC_INF for none) */
    int flags;
}
padic_radix_ctx_struct;

typedef padic_radix_ctx_struct padic_radix_ctx_t[1];

#define GR_PADIC_RADIX_CTX(ctx)  ((padic_radix_ctx_struct *) (GR_CTX_DATA_AS_PTR(ctx)))
#define PADIC_RADIX_CTX_RADIX(ctx)  (&(GR_PADIC_RADIX_CTX(ctx)->radix))
#define PADIC_RADIX_CTX_PREC_ABS(ctx) (GR_PADIC_RADIX_CTX(ctx)->prec_abs)
#define PADIC_RADIX_CTX_PREC_REL(ctx) (GR_PADIC_RADIX_CTX(ctx)->prec_rel)
#define PADIC_RADIX_CTX_FLAGS(ctx)  (GR_PADIC_RADIX_CTX(ctx)->flags)
#define PADIC_RADIX_CTX_SIGNED(ctx) (GR_PADIC_RADIX_CTX(ctx)->flags & PADIC_RADIX_SIGNED)
#define PADIC_RADIX_CTX_DECIMAL(ctx) (GR_PADIC_RADIX_CTX(ctx)->flags & PADIC_RADIX_DECIMAL)

/* Context management */
int gr_ctx_init_padic_radix(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags);
void padic_radix_ctx_clear(gr_ctx_t ctx);
int padic_radix_ctx_write(gr_stream_t out, gr_ctx_t ctx);

/* Memory management */
void padic_radix_init(padic_radix_t res, gr_ctx_t ctx);
void padic_radix_clear(padic_radix_t res, gr_ctx_t ctx);
void padic_radix_swap(padic_radix_t x, padic_radix_t y, gr_ctx_t ctx);
void padic_radix_set_shallow(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);

/* Error term access */
slong padic_radix_get_error(const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_is_exact(const padic_radix_t x, gr_ctx_t ctx);

int _padic_radix_finalize(padic_radix_t res, gr_ctx_t ctx);

/* Basic assignment and predicates */
int padic_radix_randtest(padic_radix_t res, flint_rand_t state, gr_ctx_t ctx);
int padic_radix_write(gr_stream_t out, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_zero(padic_radix_t res, gr_ctx_t ctx);
int padic_radix_one(padic_radix_t res, gr_ctx_t ctx);
int padic_radix_set(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx);
int padic_radix_set_si(padic_radix_t res, slong x, gr_ctx_t ctx);
int padic_radix_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx);
int padic_radix_exact_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx);
int padic_radix_exact_set_si(padic_radix_t res, slong x, gr_ctx_t ctx);
int padic_radix_exact_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx);
int padic_radix_get_fmpz(fmpz_t res, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_neg(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_add(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);
int padic_radix_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);
int padic_radix_mul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);
int padic_radix_inv(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_div(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);
int padic_radix_sqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);
int padic_radix_rsqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_is_square(const padic_radix_t x, gr_ctx_t ctx);

truth_t padic_radix_is_zero(const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_is_one(const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_is_neg_one(const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_is_invertible(const padic_radix_t x, gr_ctx_t ctx);
truth_t padic_radix_equal(const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);

int padic_radix_dot(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1,
    const padic_radix_struct * vec2, slong len, gr_ctx_t ctx);
int padic_radix_dot_rev(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1,
    const padic_radix_struct * vec2, slong len, gr_ctx_t ctx);

int padic_radix_dot_strided(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);
int padic_radix_dot_strided_delayed(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);
int padic_radix_dot_strided_naive(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx);

/* For testing */
int _padic_radix_add_sub_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx);
int _padic_radix_mul_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
