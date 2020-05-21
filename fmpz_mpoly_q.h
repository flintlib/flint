/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_Q_H
#define FMPZ_MPOLY_Q_H

#ifdef FMPZ_MPOLY_Q_INLINES_C
#define FMPZ_MPOLY_Q_INLINE
#else
#define FMPZ_MPOLY_Q_INLINE static __inline__
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "flint/fmpz_mpoly.h"

typedef struct
{
    fmpz_mpoly_struct num;
    fmpz_mpoly_struct den;
}
fmpz_mpoly_q_struct;

typedef fmpz_mpoly_q_struct fmpz_mpoly_q_t[1];

#define fmpz_mpoly_q_numref(x) (&((x)->num))
#define fmpz_mpoly_q_denref(x) (&((x)->den))

/* Memory management */

void fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx);

/* Assignment */

void fmpz_mpoly_q_swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx);

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

/* Polynomial helper functions */

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_gcd_assert_successful(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (!fmpz_mpoly_gcd(res, x, y, ctx))
    {
        flint_printf("fmpz_mpoly_gcd failed\n");
        flint_abort();
    }
}

FMPZ_MPOLY_Q_INLINE void
fmpz_mpoly_divexact(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(y, ctx))
        fmpz_mpoly_scalar_divexact_fmpz(res, x, y->coeffs, ctx);
    else
        fmpz_mpoly_div(res, x, y, ctx);
}

#ifdef __cplusplus
}
#endif

#endif

