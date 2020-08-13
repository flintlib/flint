/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

ulong
_fmpz_poly_is_cyclotomic(const fmpz * poly, slong len)
{
    ulong * phi;
    ulong p, q, N1, N2;
    slong i, d;
    ulong res;
    double U;
    fmpz_poly_t tmp;

    d = len - 1;

    if (d < 1)
        return 0;

    if (d == 1)
    {
        if (fmpz_is_one(poly + 1))
        {
            if (fmpz_is_one(poly))
                return 2;
            if (fmpz_equal_si(poly, -1))
                return 1;
        }

        return 0;
    }

    if (d % 2 != 0)
        return 0;

    if (!fmpz_is_one(poly))
        return 0;

    /* Rule out non-palindromes */
    for (i = 0; i < d / 2; i++)
    {
        if (!fmpz_equal(poly + i, poly + d - i))
            return 0;
    }

    /* Compute inverse image of the totient function to find candidate
       indices of the cyclotomic polynomial */

    /* Determine lower and upper bounds [N1, N2) */
    U = d;
    for (p = 2; p <= d; p++)
        if (d % (p - 1) == 0 && n_is_prime(p))
            U = (U * p) / (p - 1);
    N1 = d + 1;
    N2 = U + 3;   /* +3 as safety for possible floating-point rounding */

    res = 0;
    phi = flint_malloc(N2 * sizeof(ulong));
    fmpz_poly_init(tmp);

    for (i = 0; i < N2; i++)
        phi[i] = i;

    for (p = 2; p < N2; p++)
    {
        if (phi[p] == p)
        {
            phi[p] = p - 1;
            for (q = 2 * p; q < N2; q += p)
                phi[q] = (phi[q] / p) * (p - 1);
        }
    }

    for (i = N1; i < N2 && !res; i++)
    {
        if (phi[i] == d)
        {
            /* todo: we could avoid O(len) overhead by computing the
               factorisation of phi and calling _fmpz_poly_cyclotomic,
               checking only the deflated polynomial */
            fmpz_poly_cyclotomic(tmp, i);
            if (tmp->length == len && _fmpz_vec_equal(poly, tmp->coeffs, len))
                res = i;
        }
    }

    flint_free(phi);
    fmpz_poly_clear(tmp);

    return res;
}

ulong
fmpz_poly_is_cyclotomic(const fmpz_poly_t poly)
{
    return _fmpz_poly_is_cyclotomic(poly->coeffs, poly->length);
}
