/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly_q.h"

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
