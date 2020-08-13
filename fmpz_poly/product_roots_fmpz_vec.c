/*
    Copyright (C) 2011 Fredrik Johansson

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

void
_fmpz_poly_product_roots_fmpz_vec(fmpz * poly, const fmpz * xs, slong n)
{
    if (n == 0)
    {
        fmpz_one(poly);
    }
    else if (n < 20)
    {
        slong i, j;

        fmpz_one(poly + n);
        fmpz_neg(poly + n - 1, xs);

        for (i = 1; i < n; i++)
        {
            fmpz_mul(poly + n - i - 1, poly + n - i, xs + i);
            fmpz_neg(poly + n - i - 1, poly + n - i - 1);
            for (j = 0; j < i - 1; j++)
                fmpz_submul(poly + n - i + j, poly + n - i + j + 1, xs + i);
            fmpz_sub(poly + n - 1, poly + n - 1, xs + i);
        }
    }
    else
    {
        slong m;
        fmpz * tmp;

        m = (n + 1) / 2;

        tmp = _fmpz_vec_init(n + 2);

        _fmpz_poly_product_roots_fmpz_vec(tmp, xs, m);
        _fmpz_poly_product_roots_fmpz_vec(tmp + m + 1, xs + m, n - m);
        _fmpz_poly_mul(poly, tmp, m + 1, tmp + m + 1, n - m + 1);

        _fmpz_vec_clear(tmp, n + 2);
    }
}

void
fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly, const fmpz * xs, slong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_product_roots_fmpz_vec(poly->coeffs, xs, n);
    _fmpz_poly_set_length(poly, n + 1);
}
