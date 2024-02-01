/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly_q.h"

void
fmpz_poly_q_scalar_div_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)
{
    fmpz_t num, den;

    if (fmpz_sgn(fmpq_numref(x)) == 0)
        flint_throw(FLINT_ERROR, "Division by zero in %s\n", __func__);

    fmpz_init(num);
    fmpz_init(den);
    fmpz_set(num, fmpq_numref(x));
    fmpz_set(den, fmpq_denref(x));

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, den);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, num);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(num);
    fmpz_clear(den);
}

void
fmpz_poly_q_scalar_div_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
{
    fmpz_t y;

    if (fmpz_sgn(x) == 0)
        flint_throw(FLINT_ERROR, "Division by zero in %s\n", __func__);

    fmpz_init(y);
    fmpz_set(y, x);

    fmpz_poly_set(rop->num, op->num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, y);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}

void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
{
    fmpz_t cont, fx, gcd;

    if (FLINT_ABS(x) <= 1)
    {
        if (x == 0)
        {
            flint_throw(FLINT_ERROR, "Exception (fmpz_poly_q_scalar_div_si). Division by zero.\n");
        }
        if (x == 1)
            fmpz_poly_q_set(rop, op);
        else
            fmpz_poly_q_neg(rop, op);
        return;
    }
    if (fmpz_poly_q_is_zero(op))
    {
        fmpz_poly_q_zero(rop);
        return;
    }

    fmpz_init(cont);
    fmpz_poly_content(cont, op->num);

    if (fmpz_is_one(cont))
    {
        if (x > 0)
        {
            fmpz_poly_set(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, x);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_poly_scalar_mul_ui(rop->den, op->den, - ((ulong) x));
        }
        fmpz_clear(cont);
        return;
    }

    fmpz_init(fx);
    fmpz_init(gcd);

    fmpz_set_si(fx, x);
    fmpz_gcd(gcd, cont, fx);

    if (fmpz_is_one(gcd))
    {
        if (x > 0)
        {
            fmpz_poly_set(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, x);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_poly_scalar_mul_ui(rop->den, op->den, - ((ulong) x));
        }
    }
    else
    {
        fmpz_poly_scalar_divexact_fmpz(rop->num, op->num, gcd);
        fmpz_divexact(fx, fx, gcd);
        fmpz_poly_scalar_mul_fmpz(rop->den, op->den, fx);
        if (x < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_poly_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}

void
fmpz_poly_q_scalar_mul_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)
{
    fmpz_t num, den;

    fmpz_init(num);
    fmpz_init(den);
    fmpz_set(num, fmpq_numref(x));
    fmpz_set(den, fmpq_denref(x));

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, den);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(num);
    fmpz_clear(den);
}

void
fmpz_poly_q_scalar_mul_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
{
    fmpz_t y;

    fmpz_init(y);
    fmpz_set(y, x);

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, y);
    fmpz_poly_set(rop->den, op->den);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}

void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
{
    fmpz_t cont, fx, gcd;

    if (fmpz_poly_q_is_zero(op) || (x == 0))
    {
        fmpz_poly_q_zero(rop);
        return;
    }

    if (x == 1)
    {
        fmpz_poly_q_set(rop, op);
        return;
    }

    fmpz_init(cont);
    fmpz_poly_content(cont, op->den);

    if (fmpz_is_one(cont))
    {
        fmpz_poly_scalar_mul_si(rop->num, op->num, x);
        fmpz_poly_set(rop->den, op->den);
        fmpz_clear(cont);
        return;
    }

    fmpz_init(fx);
    fmpz_init(gcd);

    fmpz_set_si(fx, x);
    fmpz_gcd(gcd, cont, fx);

    if (fmpz_is_one(gcd))
    {
        fmpz_poly_scalar_mul_si(rop->num, op->num, x);
        fmpz_poly_set(rop->den, op->den);
    }
    else
    {
        fmpz_divexact(fx, fx, gcd);
        fmpz_poly_scalar_mul_fmpz(rop->num, op->num, fx);
        fmpz_poly_scalar_divexact_fmpz(rop->den, op->den, gcd);
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}
