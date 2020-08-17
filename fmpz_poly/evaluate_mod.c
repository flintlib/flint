/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

mp_limb_t _fmpz_poly_evaluate_mod(const fmpz * poly, slong len, mp_limb_t a, 
                                  mp_limb_t n, mp_limb_t ninv)
{
    mp_limb_t c, res = 0;

    while (len--)
    {
        c = fmpz_fdiv_ui(poly + len, n);
        res = n_addmod(n_mulmod2_preinv(res, a, n, ninv), c, n);
    }

    return res;
}

mp_limb_t fmpz_poly_evaluate_mod(const fmpz_poly_t poly, mp_limb_t a, 
                                 mp_limb_t n)
{
    if (poly->length == 0)
        return 0;

    if (a == 0)
    {
        mp_limb_t res;
        res = fmpz_fdiv_ui(poly->coeffs, n);
        return res;
    }
    else
    {
        mp_limb_t ninv;

        ninv = n_preinvert_limb(n);
        return _fmpz_poly_evaluate_mod(poly->coeffs, poly->length, a, n, ninv);
    }
}

