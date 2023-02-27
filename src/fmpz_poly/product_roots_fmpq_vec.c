/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpq.h"
#include "fmpz_poly.h"

FLINT_DLL void _fmpz_poly_product_roots_fmpq_vec(fmpz * poly,
                                        const fmpq * xs, slong n)
{
    if (n == 0)
    {
        fmpz_one(poly);
    }
    else if (n < 20)
    {
        slong i,j;

        fmpz_set(poly + n, fmpq_denref(xs));
        fmpz_set(poly + n - 1, fmpq_numref(xs));
        fmpz_neg(poly + n - 1, poly + n - 1);

        for (i = 1; i < n; i++)
        {
            fmpz_mul(poly + n - i - 1, poly + n - i, fmpq_numref(xs + i));
            fmpz_neg(poly + n - i - 1, poly + n - i - 1);
            for (j = 0; j < i; j++)
            {
                fmpz_mul(poly + n - i + j, poly + n - i + j, fmpq_denref(xs + i));
                fmpz_submul(poly + n - i + j, poly + n - i + j + 1, fmpq_numref(xs + i));
            }
            fmpz_mul(poly + n, poly + n, fmpq_denref(xs + i));
        }
    }
    else
    {
        slong m;
        fmpz * tmp;

        m = (n + 1) / 2;

        tmp = _fmpz_vec_init(n + 2);

        _fmpz_poly_product_roots_fmpq_vec(tmp, xs, m);
        _fmpz_poly_product_roots_fmpq_vec(tmp + m + 1, xs + m, n - m);
        _fmpz_poly_mul(poly, tmp, m + 1, tmp + m + 1, n - m + 1);

        _fmpz_vec_clear(tmp, n + 2);
    }
}

void fmpz_poly_product_roots_fmpq_vec(fmpz_poly_t poly,
                                        const fmpq * xs, slong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_product_roots_fmpq_vec(poly->coeffs, xs, n);
    _fmpz_poly_set_length(poly, n + 1);
}
