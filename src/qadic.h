/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QADIC_H
#define QADIC_H

#ifdef QADIC_INLINES_C
#define QADIC_INLINE
#else
#define QADIC_INLINE static __inline__
#endif

#include "padic_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Data types and context ****************************************************/

typedef padic_poly_struct qadic_struct;
typedef padic_poly_t qadic_t;

QADIC_INLINE slong qadic_val(const qadic_t op) { return padic_poly_val(op); }
QADIC_INLINE slong qadic_prec(const qadic_t op) { return padic_poly_prec(op); }

typedef struct
{
    padic_ctx_struct pctx;

    fmpz * a;
    slong * j;
    slong len;

    char * var;
}
qadic_ctx_struct;

typedef qadic_ctx_struct qadic_ctx_t[1];

void qadic_ctx_init(
    qadic_ctx_t ctx,
    const fmpz_t p,
    slong d,
    slong min,
    slong max,
    const char * var,
    enum padic_print_mode mode);
void qadic_ctx_init_conway(
    qadic_ctx_t ctx,
    const fmpz_t p,
    slong d,
    slong min,
    slong max,
    const char * var,
    enum padic_print_mode mode);
void qadic_ctx_clear(qadic_ctx_t ctx);

QADIC_INLINE
slong qadic_ctx_degree(const qadic_ctx_t ctx)
{
    return ctx->j[ctx->len - 1];
}

void qadic_ctx_print(const qadic_ctx_t ctx);

/* Memory management *********************************************************/

QADIC_INLINE void qadic_init(qadic_t x) { padic_poly_init(x); }
QADIC_INLINE
void qadic_init2(qadic_t rop, slong prec)
{
    padic_poly_init2(rop, 0, prec);
}
QADIC_INLINE void qadic_clear(qadic_t x) { padic_poly_clear(x); }

/* TODO: Consider renaming these functions, prefix for the "qadic" module. */
void _fmpz_poly_reduce(fmpz * R, slong lenR, const fmpz * a, const slong * j, slong len);
void _fmpz_mod_poly_reduce(
    fmpz * R,
    slong lenR,
    const fmpz * a,
    const slong * j,
    slong len,
    const fmpz_t p);

void qadic_reduce(qadic_t x, const qadic_ctx_t ctx);

/* Randomisation *************************************************************/

QADIC_INLINE
void qadic_randtest(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    padic_poly_randtest(x, state, qadic_ctx_degree(ctx), &ctx->pctx);
}
QADIC_INLINE
void qadic_randtest_not_zero(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    padic_poly_randtest_not_zero(x, state, qadic_ctx_degree(ctx), &ctx->pctx);
}
QADIC_INLINE
void qadic_randtest_val(qadic_t x, flint_rand_t state, slong val, const qadic_ctx_t ctx)
{
    padic_poly_randtest_val(x, state, val, qadic_ctx_degree(ctx), &ctx->pctx);
}

QADIC_INLINE
void qadic_randtest_int(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(x);

    if (N <= 0)
    {
        padic_poly_zero(x);
    }
    else
    {
        padic_poly_randtest_val(x, state, n_randint(state, N), qadic_ctx_degree(ctx), &ctx->pctx);
    }
}

/* Assignments and conversions ***********************************************/

QADIC_INLINE void qadic_zero(qadic_t op) { padic_poly_zero(op); }
QADIC_INLINE void qadic_one(qadic_t op) { padic_poly_one(op); }

void qadic_gen(qadic_t x, const qadic_ctx_t ctx);

QADIC_INLINE
void qadic_set_ui(qadic_t rop, ulong op, const qadic_ctx_t ctx)
{
    padic_poly_set_ui(rop, op, &ctx->pctx);
}
void qadic_set_fmpz_poly(qadic_t rop, const fmpz_poly_t op, const qadic_ctx_t ctx);
QADIC_INLINE
void qadic_set(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    padic_poly_set(rop, op, &(ctx->pctx));
}

int qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

/* Comparison ****************************************************************/

QADIC_INLINE int qadic_is_zero(const qadic_t op) { return padic_poly_is_zero(op); }
QADIC_INLINE int qadic_is_one(const qadic_t op) { return padic_poly_is_one(op); }

QADIC_INLINE
int qadic_equal(const qadic_t op1, const qadic_t op2)
{
    return padic_poly_equal(op1, op2);
}

/* Basic arithmetic **********************************************************/

QADIC_INLINE
void qadic_add(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
{
    padic_poly_add(x, y, z, &ctx->pctx);
}
QADIC_INLINE
void qadic_sub(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
{
    padic_poly_sub(x, y, z, &ctx->pctx);
}
QADIC_INLINE
void qadic_neg(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
{
    padic_poly_neg(x, y, &ctx->pctx);
}

void qadic_mul(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx);

void _qadic_inv(
    fmpz * rop,
    const fmpz * op,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx);

void _qadic_pow(
    fmpz * rop,
    const fmpz * op,
    slong len,
    const fmpz_t e,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p);
void qadic_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx);

/* Special functions *********************************************************/

void _qadic_exp_rectangular(
    fmpz * rop,
    const fmpz * op,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_exp_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_exp_balanced(
    fmpz * rop,
    const fmpz * op,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_exp_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_exp(
    fmpz * rop,
    const fmpz * op,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log_rectangular(
    fmpz * z,
    const fmpz * y,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log_balanced(
    fmpz * z,
    const fmpz * y,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_log_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log(
    fmpz * z,
    const fmpz * y,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N,
    const fmpz_t pN);
int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_frobenius_a(
    fmpz * rop,
    slong exp,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void _qadic_frobenius(
    fmpz * rop,
    const fmpz * op,
    slong len,
    slong e,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_frobenius(qadic_t rop, const qadic_t op, slong e, const qadic_ctx_t ctx);

void _qadic_teichmuller(
    fmpz * rop,
    const fmpz * op,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_trace(
    fmpz_t rop,
    const fmpz * op,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t pN);
void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_norm_resultant(
    fmpz_t rop,
    const fmpz * op,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_norm_analytic(
    fmpz_t rop,
    const fmpz * y,
    slong v,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_norm(
    fmpz_t rop,
    const fmpz * op,
    slong len,
    const fmpz * a,
    const slong * j,
    slong lena,
    const fmpz_t p,
    slong N);
void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

/* Output ********************************************************************/

#ifdef FLINT_HAVE_FILE
int qadic_fprint_pretty(FILE * file, const qadic_t op, const qadic_ctx_t ctx);
#endif

int qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx);
int qadic_debug(const qadic_t op);

#ifdef __cplusplus
}
#endif

#endif
