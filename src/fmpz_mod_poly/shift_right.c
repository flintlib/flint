/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009, 2008 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong i;

    /* Copy in forward order to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = 0; i < len - n; i++)
            fmpz_set(res + i, poly + n + i);
    }
    else
    {
        for (i = 0; i < len - n; i++)
            fmpz_swap(res + i, res + n + i);
    }

}

void fmpz_mod_poly_shift_right(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t ctx)
{
    if (n == 0)
    {
        fmpz_mod_poly_set(res, poly, ctx);
        return;
    }

    if (poly->length <= n)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(res, poly->length - n, ctx);
    _fmpz_mod_poly_shift_right(res->coeffs, poly->coeffs, poly->length, n);
    _fmpz_mod_poly_set_length(res, poly->length - n);
}
