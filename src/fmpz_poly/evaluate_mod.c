/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

ulong _fmpz_poly_evaluate_mod(const fmpz * poly, slong len, ulong a,
                                  ulong n, ulong ninv)
{
    ulong c, res = 0;

    while (len--)
    {
        c = fmpz_fdiv_ui(poly + len, n);
        res = n_addmod(n_mulmod2_preinv(res, a, n, ninv), c, n);
    }

    return res;
}

ulong fmpz_poly_evaluate_mod(const fmpz_poly_t poly, ulong a,
                                 ulong n)
{
    if (poly->length == 0)
        return 0;

    if (a == 0)
    {
        ulong res;
        res = fmpz_fdiv_ui(poly->coeffs, n);
        return res;
    }
    else
    {
        ulong ninv;

        ninv = n_preinvert_limb(n);
        return _fmpz_poly_evaluate_mod(poly->coeffs, poly->length, a, n, ninv);
    }
}
