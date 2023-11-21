/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

#if FLINT_BITS == 64
#define GCD_MAX_D WORD(1125899906842623)
#define GCD_MIN_D WORD(-1125899906842623)
#else
#define GCD_MAX_D COEFF_MAX
#define GCD_MIN_D COEFF_MIN
#endif

void
_fmpzi_gcd_dddd(fmpzi_t res, double a, double b, double c, double d);

void fmpzi_gcd(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpz a, b, c, d;
    slong bits1, bits2;

    if (fmpzi_is_zero(x))
    {
        fmpzi_canonicalise_unit(res, y);
        return;
    }

    if (fmpzi_is_zero(y))
    {
        fmpzi_canonicalise_unit(res, x);
        return;
    }

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
        return;
    }

    if ((!COEFF_IS_MPZ(a) && !COEFF_IS_MPZ(b)) ||
        (!COEFF_IS_MPZ(c) && !COEFF_IS_MPZ(d)))
    {
        fmpzi_gcd_euclidean_improved(res, x, y);
        return;
    }

    bits1 = fmpzi_bits(x);
    bits2 = fmpzi_bits(y);

    if (bits1 <= 30000 || bits2 <= 30000)
    {
        fmpzi_gcd_euclidean_improved(res, x, y);
        return;
    }

    if (bits1 > bits2 * 1.5 + 64)
    {
        fmpzi_t q;
        fmpzi_init(q);
        fmpzi_divrem_approx(q, res, x, y);
        fmpzi_gcd_shortest(res, res, y);
        fmpzi_clear(q);
    }
    else if (bits2 > bits1 * 1.5 + 64)
    {
        fmpzi_t q;
        fmpzi_init(q);
        fmpzi_divrem_approx(q, res, y, x);
        fmpzi_gcd_shortest(res, x, res);
        fmpzi_clear(q);
    }
    else
    {
        fmpzi_gcd_shortest(res, x, y);
    }
}
