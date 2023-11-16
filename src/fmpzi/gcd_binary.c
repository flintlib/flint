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

double
fmpzi_norm_approx_d_2exp(slong * exp, const fmpzi_t x)
{
    double a, b;
    slong aexp, bexp;
    int e;

    a = fmpz_get_d_2exp(&aexp, fmpzi_realref(x));
    b = fmpz_get_d_2exp(&bexp, fmpzi_imagref(x));

    if (aexp >= bexp)
    {
        if (aexp >= bexp + 64)
            b = 0.0;
        else
            b = ldexp(b, aexp - bexp);
    }
    else
    {
        if (bexp >= aexp + 64)
            a = 0.0;
        else
            a = ldexp(a, bexp - aexp);
    }

    a = a * a + b * b;
    a = frexp(a, &e);
    aexp += e;

    *exp = aexp;
    return a;
}

double
fmpzi_norm_approx_d(const fmpzi_t x)
{
    double a, b;

    a = fmpz_get_d(fmpzi_realref(x));
    b = fmpz_get_d(fmpzi_imagref(x));

    return a * a + b * b;
}

void
fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t X, const fmpzi_t Y)
{
    fmpzi_t x, y, z;
    slong hx, hy;

    if (fmpzi_is_zero(X))
    {
        fmpzi_canonicalise_unit(res, Y);
        return;
    }

    if (fmpzi_is_zero(Y))
    {
        fmpzi_canonicalise_unit(res, X);
        return;
    }

    /* not implemented */
    if (fmpzi_bits(X) > 500 || fmpzi_bits(Y) > 500)
    {
        fmpzi_gcd_euclidean(res, X, Y);
        return;
    }

    fmpzi_init(x);
    fmpzi_init(y);
    fmpzi_init(z);

    hx = fmpzi_remove_one_plus_i(x, X);
    hy = fmpzi_remove_one_plus_i(y, Y);

    if (fmpzi_norm_approx_d(x) < fmpzi_norm_approx_d(y))
        fmpzi_swap(x, y);

    while (!fmpzi_is_zero(y))
    {
        double a, b, c, d;
        double N, N1, N2, N3, N4;

        a = fmpz_get_d(fmpzi_realref(x));
        b = fmpz_get_d(fmpzi_imagref(x));
        c = fmpz_get_d(fmpzi_realref(y));
        d = fmpz_get_d(fmpzi_imagref(y));

        N1 = (a + b) * (a + b) + (c + d) * (c + d);
        N2 = (a - b) * (a - b) + (c - d) * (c - d);
        N3 = (a - d) * (a - d) + (b + c) * (b + c);
        N4 = (a + d) * (a + d) + (b - c) * (b - c);

        N = FLINT_MIN(FLINT_MIN(N1, N2), FLINT_MIN(N3, N4));

        if (N == N1)
        {
            fmpzi_add(z, x, y);
        }
        else if (N == N2)
        {
            fmpzi_sub(z, x, y);
        }
        else if (N == N3)
        {
            fmpz_sub(fmpzi_realref(z), fmpzi_realref(x), fmpzi_imagref(y));
            fmpz_add(fmpzi_imagref(z), fmpzi_imagref(x), fmpzi_realref(y));
        }
        else
        {
            fmpz_add(fmpzi_realref(z), fmpzi_realref(x), fmpzi_imagref(y));
            fmpz_sub(fmpzi_imagref(z), fmpzi_imagref(x), fmpzi_realref(y));
        }

        fmpzi_remove_one_plus_i(x, z);

        if (fmpzi_norm_approx_d(x) < fmpzi_norm_approx_d(y))
            fmpzi_swap(x, y);
    }

    fmpzi_swap(res, x);

    hx = FLINT_MIN(hx, hy);
    if (hx != 0)
    {
        fmpzi_set_si_si(x, 1, 1);
        fmpzi_pow_ui(x, x, hx);
        fmpzi_mul(res, res, x);
    }

    fmpzi_canonicalise_unit(res, res);

    fmpzi_clear(x);
    fmpzi_clear(y);
    fmpzi_clear(z);
}
