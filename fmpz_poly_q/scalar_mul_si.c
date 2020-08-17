/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

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
