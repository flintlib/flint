/*
   Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int _fmpz_poly_has_real_root(const fmpz * p, slong len)
{
    slong n, i;
    int s, t;

    /* O(1) conditions:                          */
    /*    - constant polynomial                  */
    /*    - odd degree                           */
    /*    - p(0) = 0                             */
    /*    - sign(p(0)) * sign(p(+infinity)) = -1 */
    if (len == 1)
        return 0;
    if (len % 2 == 0)
        return 1;
    if (fmpz_is_zero(p) || (fmpz_sgn(p) * fmpz_sgn(p + len - 1) < 0))
       return 1;

    /* O(len) conditions: Descartes rule of sign */
    n = 0;
    s = fmpz_sgn(p);
    for (i = 1; i < len; i++)
    {
        if (fmpz_is_zero(p + i)) continue;
        t = fmpz_sgn(p + i);
        if (t != s)
        {
            n += 1;
            s = t;
        }
    }
    if (n % 2 == 1) return 1;

    n = 0;
    s = fmpz_sgn(p);
    for (i = 1; i < len; i++)
    {
        if (fmpz_is_zero(p + i)) continue;
        t = fmpz_sgn(p + i);
        if (i % 2) t = -t;
        if (t != s)
        {
            n += 1;
            s = t;
        }
    }
    if (n % 2 == 1) return 1;

    /* try to isolate one root */
    return _fmpz_poly_num_real_roots(p, len) != 0;
}

int fmpz_poly_has_real_root(const fmpz_poly_t pol)
{
    return _fmpz_poly_has_real_root(pol->coeffs, pol->length);
}

