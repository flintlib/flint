/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e)
{
    ulong i, f;
    fmpz_t a, b, c;

    *a = WORD(1);
    *b = WORD(1);
    *c = WORD(1);

    fmpz_one(res);
    fmpz_one(res + e);

    for (i = UWORD(1), f = e - UWORD(1); i <= (e - UWORD(1)) >> 1; i++, f--)
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);
        fmpz_mul_ui(c, c, f + UWORD(1));
        fmpz_divexact_ui(c, c, i);

        fmpz_mul(res + i, b, c);
        fmpz_mul(res + f, a, c);
    }

    if ((e & UWORD(1)) == UWORD(0))
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);
        fmpz_mul_ui(c, c, f + UWORD(1));
        fmpz_divexact_ui(c, c, i);

        fmpz_mul(res + i, b, c);
        fmpz_mul(res + i, res + i, a);
        i++, f--;
    }

    for ( ; i <= e; i++, f--)
    {
        fmpz_mul(a, a, poly);
        fmpz_mul(b, b, poly + 1);

        fmpz_mul(res + i, res + i, b);
        fmpz_mul(res + f, res + f, a);
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

void
fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const slong len = poly->length;
    slong rlen;

    if (len != 2)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_pow_binomial). poly->length not equal to 2.\n");
    }

    if (e < UWORD(3))
    {
        if (e == UWORD(0))
            fmpz_poly_set_ui(res, UWORD(1));
        else if (e == UWORD(1))
            fmpz_poly_set(res, poly);
        else  /* e == UWORD(2) */
            fmpz_poly_sqr(res, poly);
        return;
    }

    rlen = (slong) e + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_set_length(res, rlen);
        _fmpz_poly_pow_binomial(res->coeffs, poly->coeffs, e);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_set_length(t, rlen);
        _fmpz_poly_pow_binomial(t->coeffs, poly->coeffs, e);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
