/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

static void
_interpolate_newton(arb_ptr ys, arb_srcptr xs, slong n, slong prec)
{
    arb_t p, q, t;
    slong i, j;

    arb_init(p);
    arb_init(q);
    arb_init(t);

    for (i = 1; i < n; i++)
    {
        arb_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            arb_sub(p, ys + j, t, prec);
            arb_sub(q, xs + j, xs + j - i, prec);
            arb_set(t, ys + j);
            arb_div(ys + j, p, q, prec);
        }
    }

    arb_clear(p);
    arb_clear(q);
    arb_clear(t);
}

static void
_newton_to_monomial(arb_ptr ys, arb_srcptr xs, slong n, slong prec)
{
    arb_t t, u;
    slong i, j;

    arb_init(t);
    arb_init(u);

    for (i = n - 2; i >= 0; i--)
    {
        arb_set(t, ys + i);
        arb_set(ys + i, ys + i + 1);

        for (j = i + 1; j < n - 1; j++)
        {
            arb_mul(u, ys + j, xs + i, prec);
            arb_sub(ys + j, ys + j + 1, u, prec);
        }

        arb_mul(u, ys + n - 1, xs + i, prec);
        arb_sub(ys + n - 1, t, u, prec);
    }

    _arb_poly_reverse(ys, ys, n, n);

    arb_clear(t);
    arb_clear(u);
}

void
_arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs,
    arb_srcptr ys, slong n, slong prec)
{
    if (n == 1)
    {
        arb_set(poly, ys);
    }
    else
    {
        _arb_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, prec);
        while (n > 0 && arb_is_zero(poly + n - 1)) n--;
        _newton_to_monomial(poly, xs, n, prec);
    }
}

void
arb_poly_interpolate_newton(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(poly);
    }
    else
    {
        arb_poly_fit_length(poly, n);
        _arb_poly_set_length(poly, n);
        _arb_poly_interpolate_newton(poly->coeffs,
            xs, ys, n, prec);
        _arb_poly_normalise(poly);
    }
}
