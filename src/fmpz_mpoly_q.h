/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_Q_H
#define FMPZ_MPOLY_Q_H

#ifdef FMPZ_MPOLY_Q_INLINES_C
#define FMPZ_MPOLY_Q_INLINE
#else
#define FMPZ_MPOLY_Q_INLINE static inline
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "acb_types.h"

#define fmpz_mpoly_q_numref(x) (&((x)->num))
#define fmpz_mpoly_q_denref(x) (&((x)->den))

/* Memory management */

void fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx);

/* Assignment */

void fmpz_mpoly_q_swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong x, const fmpz_mpoly_ctx_t ctx);

/* Canonicalisation */

void fmpz_mpoly_q_canonicalise(fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx);

int fmpz_mpoly_q_is_canonical(const fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx);

/* Properties */

FMPZ_MPOLY_Q_INLINE int
fmpz_mpoly_q_is_zero(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), ctx);
}

FMPZ_MPOLY_Q_INLINE int
fmpz_mpoly_q_is_one(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_one(fmpz_mpoly_q_numref(x), ctx) &&
           fmpz_mpoly_is_one(fmpz_mpoly_q_denref(x), ctx);
}

FMPZ_MPOLY_Q_INLINE int
fmpz_mpoly_q_is_fmpz(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(x), ctx) &&
           fmpz_mpoly_is_one(fmpz_mpoly_q_denref(x), ctx);
}

FMPZ_MPOLY_Q_INLINE int
fmpz_mpoly_q_is_fmpq(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(x), ctx) &&
           fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), ctx);
}

void fmpz_mpoly_q_used_vars(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_used_vars_num(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_used_vars_den(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx);

/* Special values */

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_zero(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_zero(fmpz_mpoly_q_numref(res), ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_one(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_one(fmpz_mpoly_q_numref(res), ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_gen(fmpz_mpoly_q_t res, slong i, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_gen(fmpz_mpoly_q_numref(res), i, ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

/* Input and output */

void fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx);

/* Random generation */

void fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, slong length, mp_limb_t coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx);

/* Comparisons */

int fmpz_mpoly_q_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);

/* Arithmetic */

void fmpz_mpoly_q_neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_add(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_sub(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_mul(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_div(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_add_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void
_fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_sub_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_mul_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_add_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mpoly_q_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_sub_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mpoly_q_sub_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_mul_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mpoly_q_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_div_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mpoly_q_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

/* Polynomial helper functions */

FMPZ_MPOLY_Q_INLINE void
_fmpz_vec_content2(fmpz_t res, const fmpz * vec, slong len, const fmpz_t inp)
{
    if (fmpz_is_pm1(inp))
    {
        fmpz_one(res);
    }
    else
    {
        slong i;
        fmpz_abs(res, inp);
        for (i = len - 1; i >= 0; i--)
        {
            fmpz_gcd(res, res, vec + i);
            if (fmpz_is_one(res))
                break;
        }
    }
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_gcd_assert_successful(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (!fmpz_mpoly_gcd(res, x, y, ctx))
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_gcd failed\n");
    }
}

FMPZ_MPOLY_Q_INLINE void
_fmpz_mpoly_q_mpoly_divexact(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(y, ctx))
        fmpz_mpoly_scalar_divexact_fmpz(res, x, y->coeffs, ctx);
    else
        fmpz_mpoly_div(res, x, y, ctx);
}

/* Content */

FMPZ_MPOLY_Q_INLINE void
_fmpz_mpoly_q_content(fmpz_t num, fmpz_t den, const fmpz_mpoly_t xnum, const fmpz_mpoly_t xden, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(xnum, ctx))
    {
        fmpz_one(num);
        fmpz_one(den);
    }
    else
    {
        _fmpz_vec_content(den, xden->coeffs, xden->length);
        _fmpz_vec_content(num, xnum->coeffs, xnum->length);
    }
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_q_content(fmpq_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_content(fmpq_numref(res), fmpq_denref(res), fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x), ctx);
}

/* Evaluation */

void fmpz_mpoly_q_evaluate_acb(acb_t res, const fmpz_mpoly_q_t f, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

