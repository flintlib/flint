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
_fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = len; i--; )
            fmpz_set(res + n + i, poly + i);
    }
    else
    {
        for (i = len; i--; )
            fmpz_swap(res + n + i, res + i);
    }

    for (i = 0; i < n; i++)
        fmpz_zero(res + i);
}

void fmpz_mod_poly_shift_left(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t ctx)
{
    if (n == 0)
    {
        fmpz_mod_poly_set(res, poly, ctx);
        return;
    }

    if (poly->length == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(res, poly->length + n, ctx);
    _fmpz_mod_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n);
    _fmpz_mod_poly_set_length(res, poly->length + n);
}

