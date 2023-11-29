/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpzi.h"

void
fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)
{
    slong xbits, ybits;

    xbits = fmpzi_bits(x);
    ybits = fmpzi_bits(y);

    if (ybits == 0)
    {
        flint_throw(FLINT_ERROR, "fmpzi_divrem_approx: division by zero\n");
    }

    if (xbits == 0)
    {
        fmpzi_zero(q);
        if (r != NULL)
            fmpzi_zero(r);
        return;
    }

    if (xbits < ybits - 2)
    {
        if (r != NULL)
            fmpzi_set(r, x);
        fmpzi_zero(q);
        return;
    }

    if (xbits < ybits + 45)
    {
        double a, b, c, d, t, u, v, w, qa, qb;
        slong aexp, bexp, cexp, dexp;

        if (xbits < 500)
        {
            a = fmpz_get_d(fmpzi_realref(x));
            b = fmpz_get_d(fmpzi_imagref(x));
            c = fmpz_get_d(fmpzi_realref(y));
            d = fmpz_get_d(fmpzi_imagref(y));
        }
        else
        {
            a = fmpz_get_d_2exp(&aexp, fmpzi_realref(x));
            b = fmpz_get_d_2exp(&bexp, fmpzi_imagref(x));
            c = fmpz_get_d_2exp(&cexp, fmpzi_realref(y));
            d = fmpz_get_d_2exp(&dexp, fmpzi_imagref(y));

            a = ldexp(a, FLINT_MAX(aexp - xbits, -1024));
            b = ldexp(b, FLINT_MAX(bexp - xbits, -1024));
            c = ldexp(c, FLINT_MAX(cexp - xbits, -1024));
            d = ldexp(d, FLINT_MAX(dexp - xbits, -1024));
        }

        t = a * c + b * d;
        u = b * c - a * d;
        v = c * c + d * d;

        w = 0.5 / v;

        t = (2.0 * t + v) * w;
        u = (2.0 * u + v) * w;

        qa = floor(t);
        qb = floor(u);

        if (r != NULL)
        {
            fmpzi_t t, u;

            fmpzi_init(t);
            fmpzi_init(u);

            fmpz_set_d(fmpzi_realref(u), qa);
            fmpz_set_d(fmpzi_imagref(u), qb);

            fmpzi_mul(t, u, y);
            fmpzi_sub(r, x, t);
            fmpzi_swap(q, u);

            fmpzi_clear(t);
            fmpzi_clear(u);
        }
        else
        {
            fmpz_set_d(fmpzi_realref(q), qa);
            fmpz_set_d(fmpzi_imagref(q), qb);
        }

        return;
    }

    fmpzi_divrem(q, r, x, y);
}
