/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart

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

void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, slong len1, 
                                      const fmpz *poly2, slong len2, 
                                      const fmpz_t p, slong n)
{
    _fmpz_poly_mullow(res, poly1, len1, poly2, len2, n);
    _fmpz_vec_scalar_mod_fmpz(res, res, n, p);
}

void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if ((len1 == 0) || (len2 == 0) || (n == 0))
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);

    if ((res == poly1) || (res == poly2))
    {
        fmpz *t = _fmpz_vec_init(n);

        if (len1 >= len2)
            _fmpz_mod_poly_mullow(t, poly1->coeffs, len1, 
                            poly2->coeffs, len2, fmpz_mod_ctx_modulus(ctx), n);
        else
            _fmpz_mod_poly_mullow(t, poly2->coeffs, len2, 
                            poly1->coeffs, len1, fmpz_mod_ctx_modulus(ctx), n);

        _fmpz_vec_clear(res->coeffs, res->alloc);
        res->coeffs = t;
        res->alloc  = n;
        res->length = n;
        _fmpz_mod_poly_normalise(res);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, n, ctx);

        if (len1 >= len2)
            _fmpz_mod_poly_mullow(res->coeffs, poly1->coeffs, len1, 
                            poly2->coeffs, len2, fmpz_mod_ctx_modulus(ctx), n);
        else
            _fmpz_mod_poly_mullow(res->coeffs, poly2->coeffs, len2, 
                            poly1->coeffs, len1, fmpz_mod_ctx_modulus(ctx), n);

        _fmpz_mod_poly_set_length(res, n);
        _fmpz_mod_poly_normalise(res);
    }
}

