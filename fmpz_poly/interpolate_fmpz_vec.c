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


static void
_fmpz_poly_interpolate_newton(fmpz * ys, const fmpz * xs, slong n)
{
    fmpz_t p, q, t, r;
    slong i, j;

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);
	fmpz_init(r);

    for (i = 1; i < n; i++)
    {
        fmpz_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            fmpz_sub(p, ys + j, t);
            fmpz_sub(q, xs + j, xs + j - i);
            fmpz_set(t, ys + j);
            fmpz_fdiv_qr(ys + j, r, p, q);
			
            if (!fmpz_is_zero(r))
            {
                fmpz_clear(r);
                fmpz_clear(t);
                fmpz_clear(q);
                fmpz_clear(p);
				
                flint_throw(FLINT_INEXACT, "Not an exact division in"     
                    "fmpz_poly_interpolate_newton");
            }
        }
    }

    fmpz_clear(r);
	fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);
}

void
fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    if (n == 0)
    {
        fmpz_poly_zero(poly);
        return;
    }
    else if (n == 1)
    {
        fmpz_poly_set_fmpz(poly, ys);
        return;
    }
    else
    {
        fmpz_poly_fit_length(poly, n);
        _fmpz_vec_set(poly->coeffs, ys, n);
        _fmpz_poly_interpolate_newton(poly->coeffs, xs, n);
        _fmpz_poly_set_length(poly, n);
        _fmpz_poly_normalise(poly);
        _fmpz_poly_newton_to_monomial(poly->coeffs, xs, poly->length);
    }
}
