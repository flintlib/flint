/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

/* Multiplication and division by scalars only touch the numerator and
   the integer constant of the denominator. */

void
fmpz_mpoly_q_factored_mul_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * num;
    fmpz_mpoly_factor_struct * den;
    fmpz_t g;

    if (fmpz_is_zero(c) || fmpz_mpoly_q_factored_is_zero(x, ctx))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
        return;
    }

    fmpz_mpoly_q_factored_set(res, x, ctx);
    num = fmpz_mpoly_q_factored_numref(res);
    den = fmpz_mpoly_q_factored_denref(res);

    if (fmpz_is_one(den->constant))
    {
        fmpz_mpoly_scalar_mul_fmpz(num, num, c, ctx);
        return;
    }

    fmpz_init(g);
    fmpz_gcd(g, c, den->constant);

    if (fmpz_is_one(g))
    {
        fmpz_mpoly_scalar_mul_fmpz(num, num, c, ctx);
    }
    else
    {
        fmpz_divexact(den->constant, den->constant, g);
        fmpz_divexact(g, c, g);
        fmpz_mpoly_scalar_mul_fmpz(num, num, g, ctx);
    }

    fmpz_clear(g);
}

void
fmpz_mpoly_q_factored_div_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * num;
    fmpz_mpoly_factor_struct * den;
    fmpz_t g, t;

    if (fmpz_is_zero(c))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_div_fmpz: division by zero\n");

    if (fmpz_mpoly_q_factored_is_zero(x, ctx))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
        return;
    }

    fmpz_mpoly_q_factored_set(res, x, ctx);
    num = fmpz_mpoly_q_factored_numref(res);
    den = fmpz_mpoly_q_factored_denref(res);

    if (fmpz_is_one(c))
        return;

    fmpz_init(g);
    fmpz_init(t);

    _fmpz_vec_content2(g, num->coeffs, num->length, c);

    if (fmpz_is_one(g))
    {
        fmpz_abs(t, c);
    }
    else
    {
        fmpz_mpoly_scalar_divexact_fmpz(num, num, g, ctx);
        fmpz_divexact(t, c, g);
        fmpz_abs(t, t);
    }

    fmpz_mul(den->constant, den->constant, t);

    if (fmpz_sgn(c) < 0)
        fmpz_mpoly_neg(num, num, ctx);

    fmpz_clear(g);
    fmpz_clear(t);
}

void
fmpz_mpoly_q_factored_mul_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpq_is_zero(c))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
        return;
    }

    fmpz_mpoly_q_factored_mul_fmpz(res, x, fmpq_numref(c), ctx);
    fmpz_mpoly_q_factored_div_fmpz(res, res, fmpq_denref(c), ctx);
}

void
fmpz_mpoly_q_factored_div_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpq_is_zero(c))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_div_fmpq: division by zero\n");

    fmpz_mpoly_q_factored_mul_fmpz(res, x, fmpq_denref(c), ctx);
    fmpz_mpoly_q_factored_div_fmpz(res, res, fmpq_numref(c), ctx);
}

/* Addition and subtraction: convert the scalar and use the general code. */

void
fmpz_mpoly_q_factored_add_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz(t, c, ctx);
    fmpz_mpoly_q_factored_add(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_sub_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz(t, c, ctx);
    fmpz_mpoly_q_factored_sub(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_add_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpq(t, c, ctx);
    fmpz_mpoly_q_factored_add(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_sub_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpq(t, c, ctx);
    fmpz_mpoly_q_factored_sub(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

/* Operations with polynomials. The polynomial is copied, since it may
   alias the numerator of res. */

void
fmpz_mpoly_q_factored_add_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly(t, c, ctx);
    fmpz_mpoly_q_factored_add(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_sub_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly(t, c, ctx);
    fmpz_mpoly_q_factored_sub(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_mul_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;
    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly(t, c, ctx);
    fmpz_mpoly_q_factored_mul(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}

void
fmpz_mpoly_q_factored_div_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_t t;

    if (fmpz_mpoly_is_zero(c, ctx))
        flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_div_fmpz_mpoly: division by zero\n");

    fmpz_mpoly_q_factored_init(t, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly(t, c, ctx);
    fmpz_mpoly_q_factored_div(res, x, t, ctx);
    fmpz_mpoly_q_factored_clear(t, ctx);
}
