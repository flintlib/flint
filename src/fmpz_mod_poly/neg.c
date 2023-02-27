/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_t p)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!fmpz_is_zero(poly + i))
            fmpz_sub(res + i, p, poly + i);
        else
            fmpz_zero(res + i);
    }
}

void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = poly->length;

    fmpz_mod_poly_fit_length(res, len, ctx);
    _fmpz_mod_poly_set_length(res, len);

    _fmpz_mod_poly_neg(res->coeffs, poly->coeffs, poly->length,
                                                    fmpz_mod_ctx_modulus(ctx));
}

