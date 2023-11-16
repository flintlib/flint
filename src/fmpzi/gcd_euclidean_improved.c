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

#if FLINT_BITS == 64
#define GCD_MAX_D WORD(1125899906842623)
#define GCD_MIN_D WORD(-1125899906842623)
#else
#define GCD_MAX_D COEFF_MAX
#define GCD_MIN_D COEFF_MIN
#endif

void
_fmpzi_gcd_dddd(fmpzi_t res, double a, double b, double c, double d)
{
    double t, u, v, w, qa, qb;

    while (c != 0 || d != 0)
    {
        t = a * c + b * d;
        u = b * c - a * d;
        v = c * c + d * d;

        w = 0.5 / v;

        qa = floor((2.0 * t + v) * w);
        qb = floor((2.0 * u + v) * w);

        t = a - (qa * c - qb * d);
        u = b - (qb * c + qa * d);

        a = c;
        b = d;

        c = t;
        d = u;
    }

    if (a < 0)
    {
        a = -a;
        b = -b;
    }

    if (b > 0 && b > a)
    {
        t = b;
        b = -a;
        a = t;
    }
    else if (b < 0 && b <= -a)
    {
        t = b;
        b = a;
        a = -t;
    }

    fmpzi_set_si_si(res, a, b);
}

void
fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t X, const fmpzi_t Y)
{
    fmpzi_t x, y, q, r;
    fmpz a, b, c, d;

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

    a = *fmpzi_realref(X);
    b = *fmpzi_imagref(X);
    c = *fmpzi_realref(Y);
    d = *fmpzi_imagref(Y);

    if (GCD_MIN_D <= a && a <= GCD_MAX_D &&
        GCD_MIN_D <= b && b <= GCD_MAX_D &&
        GCD_MIN_D <= c && c <= GCD_MAX_D &&
        GCD_MIN_D <= d && d <= GCD_MAX_D)
    {
        _fmpzi_gcd_dddd(res, a, b, c, d);
        return;
    }

    fmpzi_init(x);
    fmpzi_init(y);
    fmpzi_init(q);
    fmpzi_init(r);

    fmpzi_set(x, X);
    fmpzi_set(y, Y);

    while (!fmpzi_is_zero(y))
    {
        a = *fmpzi_realref(x);
        b = *fmpzi_imagref(x);
        c = *fmpzi_realref(y);
        d = *fmpzi_imagref(y);

        if (GCD_MIN_D <= a && a <= GCD_MAX_D &&
            GCD_MIN_D <= b && b <= GCD_MAX_D &&
            GCD_MIN_D <= c && c <= GCD_MAX_D &&
            GCD_MIN_D <= d && d <= GCD_MAX_D)
        {
            _fmpzi_gcd_dddd(res, a, b, c, d);
            goto cleanup;
        }

        fmpzi_divrem_approx(q, r, x, y);
        fmpzi_swap(x, y);
        fmpzi_swap(y, r);
    }

    fmpzi_swap(res, x);
    fmpzi_canonicalise_unit(res, res);

cleanup:
    fmpzi_clear(x);
    fmpzi_clear(y);
    fmpzi_clear(q);
    fmpzi_clear(r);
}
