/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_log1p_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
{
    arb_t a;

    flen = FLINT_MIN(flen, n);

    arb_init(a);
    arb_log1p(a, f, prec);

    if (flen == 1)
    {
        _arb_vec_zero(res + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_add_ui(res, f + 0, 1, prec);
        arb_div(res + 1, f + 1, res + 0, prec);
    }
    else if (_arb_vec_is_zero(f + 1, flen - 2))  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        arb_add_ui(res, f + 0, 1, prec);

        for (i = 1, j = d; j < n; j += d, i++)
        {
            if (i == 1)
                arb_div(res + j, f + d, res, prec);
            else
                arb_mul(res + j, res + j - d, res + d, prec);
            _arb_vec_zero(res + j - d + 1, flen - 2);
        }
        _arb_vec_zero(res + j - d + 1, n - (j - d + 1));

        for (i = 2, j = 2 * d; j < n; j += d, i++)
            arb_div_si(res + j, res + j, i % 2 ? i : -i, prec);
    }
    else
    {
        arb_ptr f_diff, f_inv;
        slong alloc;

        alloc = n + flen;
        f_inv = _arb_vec_init(alloc);
        f_diff = f_inv + n;

        arb_add_ui(f_diff, f, 1, prec);
        _arb_vec_set(f_diff + 1, f + 1, flen - 1);
        _arb_poly_inv_series(f_inv, f_diff, flen, n, prec);
        _arb_poly_derivative(f_diff, f, flen, prec);
        _arb_poly_mullow(res, f_inv, n - 1, f_diff, flen - 1, n - 1, prec);
        _arb_poly_integral(res, res, n, prec);

        _arb_vec_clear(f_inv, alloc);
    }

    arb_swap(res, a);
    arb_clear(a);
}

void
arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    slong flen = f->length;

    if (flen == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (flen == 1 /*&& !arb_contains_si(f->coeffs, -1)*/)
        n = 1;

    arb_poly_fit_length(res, n);
    _arb_poly_log1p_series(res->coeffs, f->coeffs, flen, n, prec);
    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

