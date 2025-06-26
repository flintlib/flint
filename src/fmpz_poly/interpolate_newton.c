/*
    Copyright (C) 2011, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int
_fmpz_poly_interpolate_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz_t p, q, t, r;
    slong i, j;
    int ok = 1;

    if (poly != ys)
        _fmpz_vec_set(poly, ys, n);

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);
	fmpz_init(r);

    for (i = 1; i < n; i++)
    {
        fmpz_set(t, poly + i - 1);

        for (j = i; j < n; j++)
        {
            fmpz_sub(p, poly + j, t);
            fmpz_sub(q, xs + j, xs + j - i);

            if (fmpz_is_zero(q))
            {
                ok = 0;
                goto cleanup;
            }

            fmpz_swap(t, poly + j);
            fmpz_fdiv_qr(poly + j, r, p, q);

            if (!fmpz_is_zero(r))
            {
                ok = 0;
                goto cleanup;
            }
        }
    }

    FMPZ_VEC_NORM(poly, n);
    _fmpz_poly_newton_to_monomial(poly, xs, n);

cleanup:
    fmpz_clear(r);
	fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);

    return ok;
}

int
fmpz_poly_interpolate_newton(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    int ok;
    fmpz_poly_fit_length(poly, n);
    ok = _fmpz_poly_interpolate_newton(poly->coeffs, xs, ys, n);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
    return ok;
}

void
_fmpz_poly_interpolate_exact_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz_t p, q, t;
    slong i, j;

    if (poly != ys)
        _fmpz_vec_set(poly, ys, n);

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);

    for (i = 1; i < n; i++)
    {
        fmpz_set(t, poly + i - 1);

        for (j = i; j < n; j++)
        {
            fmpz_sub(p, poly + j, t);
            fmpz_sub(q, xs + j, xs + j - i);
            fmpz_swap(t, poly + j);
            fmpz_divexact(poly + j, p, q);
        }
    }

    FMPZ_VEC_NORM(poly, n);
    _fmpz_poly_newton_to_monomial(poly, xs, n);

	fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);
}

void
fmpz_poly_interpolate_exact_newton(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz_poly_fit_length(poly, n);
    _fmpz_poly_interpolate_exact_newton(poly->coeffs, xs, ys, n);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
}

