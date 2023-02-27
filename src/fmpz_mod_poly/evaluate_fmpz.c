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
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, slong len, 
                                  const fmpz_t a, const fmpz_t p)
{
    if (len == 0)
    {
        fmpz_zero(res);
    }
    else if (len == 1 || fmpz_is_zero(a))
    {
        fmpz_set(res, poly);
    }
    else
    {
        slong i = len - 1;
        fmpz_t t;

        fmpz_init(t);
        fmpz_set(res, poly + i);
        for (i = len - 2; i >= 0; i--)
        {
            fmpz_mul(t, res, a);
            fmpz_mod(t, t, p);
            fmpz_add(res, poly + i, t);
        }
        fmpz_clear(t);

        if (fmpz_cmpabs(res, p) >= 0)
            fmpz_sub(res, res, p);
    }
}

void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz_mod_poly_t poly,
                                      const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    if (res == a)
    {
        fmpz_t t;

        fmpz_init(t);
        _fmpz_mod_poly_evaluate_fmpz(t, poly->coeffs, poly->length, 
                                     a, fmpz_mod_ctx_modulus(ctx));
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
    {
        _fmpz_mod_poly_evaluate_fmpz(res, poly->coeffs, poly->length, 
                                     a, fmpz_mod_ctx_modulus(ctx));
    }
}

