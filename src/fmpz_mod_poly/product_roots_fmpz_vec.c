/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_product_roots_fmpz_vec(fmpz * poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    if (n == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
            fmpz_zero(poly);
        else
            fmpz_one(poly);
    }
    else if (n < 20)
    {
        slong i, j;

        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
            fmpz_zero(poly + n);
        else
            fmpz_one(poly + n);

        fmpz_mod_neg(poly + n - 1, xs, ctx);

        for (i = 1; i < n; i++)
        {
            fmpz_mod_mul(poly + n - i - 1, poly + n - i, xs + i, ctx);
            fmpz_mod_neg(poly + n - i - 1, poly + n - i - 1, ctx);

            for (j = 0; j < i - 1; j++)
            {
                fmpz_submul(poly + n - i + j, poly + n - i + j + 1, xs + i);
                fmpz_mod_set_fmpz(poly + n - i + j, poly + n - i + j, ctx);
            }
            fmpz_mod_sub(poly + n - 1, poly + n - 1, xs + i, ctx);
        }
    }
    else
    {
        slong m;
        fmpz * tmp;

        m = (n + 1) / 2;

        tmp = _fmpz_vec_init(n + 2);

        _fmpz_mod_poly_product_roots_fmpz_vec(tmp, xs, m, ctx);
        _fmpz_mod_poly_product_roots_fmpz_vec(tmp + m + 1, xs + m, n - m, ctx);
        _fmpz_mod_poly_mul(poly, tmp, m + 1, tmp + m + 1, n - m + 1, ctx);

        _fmpz_vec_clear(tmp, n + 2);
    }
}

void fmpz_mod_poly_product_roots_fmpz_vec(
   fmpz_mod_poly_t poly,
   const fmpz * xs, slong xlen,
   const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, xlen + 1, ctx);
    _fmpz_mod_poly_product_roots_fmpz_vec(poly->coeffs, xs, xlen, ctx);
    _fmpz_mod_poly_set_length(poly, xlen + 1);
}
