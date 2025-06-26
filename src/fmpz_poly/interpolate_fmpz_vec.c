/*
    Copyright (C) 2011, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int
_fmpz_poly_interpolate(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
{
    slong estimate_bits;
    slong xbits, ybits;

    if (n == 0)
        return 1;

    if (n == 1)
    {
        fmpz_set(poly, ys);
        return 1;
    }

    if (n <= 25)
    {
        return _fmpz_poly_interpolate_newton(poly, xs, ys, n);
    }

    xbits = _fmpz_vec_max_bits(xs, n);
    ybits = _fmpz_vec_max_bits(ys, n);

    xbits = FLINT_ABS(xbits);
    ybits = FLINT_ABS(ybits);

    /* Guess that max_i |y_i| ~= ((max_i) |x_i|)^(n-1) * max_i |f_i| */
    estimate_bits = ybits - (n - 1) * FLINT_MAX(xbits - 1, 0);
    estimate_bits = FLINT_MAX(estimate_bits, 0);

    if (estimate_bits < FLINT_BITS - 2 || n > 10.0 * sqrt(estimate_bits))
        return _fmpz_poly_interpolate_multi_mod(poly, xs, ys, n);
    else
        return _fmpz_poly_interpolate_newton(poly, xs, ys, n);
}

int
fmpz_poly_interpolate(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    int ok;
    fmpz_poly_fit_length(poly, n);
    ok = _fmpz_poly_interpolate(poly->coeffs, xs, ys, n);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
    return ok;
}

void
_fmpz_poly_interpolate_exact(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
{
    slong estimate_bits;
    slong xbits, ybits;

    if (n == 0)
        return;

    if (n == 1)
    {
        fmpz_set(poly, ys);
        return;
    }

    if (n <= 70)
    {
        _fmpz_poly_interpolate_exact_newton(poly, xs, ys, n);
        return;
    }

    xbits = _fmpz_vec_max_bits(xs, n);
    ybits = _fmpz_vec_max_bits(ys, n);

    xbits = FLINT_ABS(xbits);
    ybits = FLINT_ABS(ybits);

    /* Guess that max_i |y_i| ~= ((max_i) |x_i|)^(n-1) * max_i |f_i| */
    estimate_bits = ybits - (n - 1) * FLINT_MAX(xbits - 1, 0);
    estimate_bits = FLINT_MAX(estimate_bits, 0);

    if (estimate_bits < FLINT_BITS - 2 || n > 18.0 * sqrt(estimate_bits))
        _fmpz_poly_interpolate_multi_mod(poly, xs, ys, n);
    else
        _fmpz_poly_interpolate_exact_newton(poly, xs, ys, n);
}

void
fmpz_poly_interpolate_exact(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz_poly_fit_length(poly, n);
    _fmpz_poly_interpolate_exact(poly->coeffs, xs, ys, n);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
}

void
fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    if (!fmpz_poly_interpolate(poly, xs, ys, n))
        flint_throw(FLINT_INEXACT, "fmpz_poly_interpolate_fmpz_vec: points not disjoint or non-integral interpolant");
}

