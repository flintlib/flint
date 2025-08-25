/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_Q_H
#define FMPZ_MOD_MPOLY_Q_H

#ifdef FMPZ_MOD_MPOLY_Q_INLINES_C
#define FMPZ_MOD_MPOLY_Q_INLINE
#else
#define FMPZ_MOD_MPOLY_Q_INLINE static inline
#endif

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mpoly.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define fmpz_mod_mpoly_q_numref(x) (&((x)->num))
#define fmpz_mod_mpoly_q_denref(x) (&((x)->den))

/* Memory management */

void fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

/* Assignment */

void fmpz_mod_mpoly_q_swap(fmpz_mod_mpoly_q_t x, fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_set_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_t x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_set_si(fmpz_mod_mpoly_q_t res, slong x, const fmpz_mod_mpoly_ctx_t ctx);

/* Canonicalisation */

void fmpz_mod_mpoly_q_canonicalise(fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_q_is_canonical(const fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

/* Properties */

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_zero(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_zero(fmpz_mod_mpoly_q_numref(x), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_one(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_numref(x), ctx) &&
           fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(x), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_fmpz_mod(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_numref(x), ctx) &&
           fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(x), ctx);
}

/* Special values */

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_zero(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_zero(fmpz_mod_mpoly_q_numref(res), ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_one(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_numref(res), ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_gen(fmpz_mod_mpoly_q_t res, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_gen(fmpz_mod_mpoly_q_numref(res), i, ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

/* Input and output */

void fmpz_mod_mpoly_q_print_pretty(const fmpz_mod_mpoly_q_t f, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);
char * fmpz_mod_mpoly_q_get_str_pretty(const fmpz_mod_mpoly_q_t f, const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_set_str_pretty(fmpz_mod_mpoly_q_t res, const char * s, const char ** vars, fmpz_mod_mpoly_ctx_t ctx);

/* Random generation */

void fmpz_mod_mpoly_q_randtest(fmpz_mod_mpoly_q_t res, flint_rand_t state, slong length, slong exp_bound, const fmpz_mod_mpoly_ctx_t ctx);

/* Comparisons */

int fmpz_mod_mpoly_q_equal(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);

/* Arithmetic */

void fmpz_mod_mpoly_q_neg(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_inv(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx);

void
_fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx);

void
_fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_t y,
            const fmpz_mod_mpoly_ctx_t ctx);

void
_fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx);

void
_fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx);

void
_fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_add_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_sub_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_mul_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_mul_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_add_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mod_mpoly_q_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_sub_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mod_mpoly_q_sub_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_mul_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, c);
    fmpz_mod_mpoly_q_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_div_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    int invertible;
    fmpz_init_set_si(t, c);
    invertible = fmpz_mod_mpoly_q_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return invertible;
}

/* Polynomial helper functions */

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_gcd_assert_successful(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t x, const fmpz_mod_mpoly_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (!fmpz_mod_mpoly_gcd(res, x, y, ctx))
    {
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_gcd failed\n");
    }
}

FMPZ_MOD_MPOLY_Q_INLINE void
_fmpz_mod_mpoly_q_mpoly_divexact(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t x, const fmpz_mod_mpoly_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t g;
    fmpz_init(g);
    
    if (fmpz_mod_mpoly_is_fmpz(y, ctx))
    {
        fmpz_mod_inv(g, y->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res, x, g, ctx);
    }
    else 
        fmpz_mod_mpoly_div(res, x, y, ctx);
    fmpz_clear(g);
}

#ifdef __cplusplus
}
#endif

#endif
