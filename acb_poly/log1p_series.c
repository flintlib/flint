/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_log1p_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
{
    acb_t a;

    flen = FLINT_MIN(flen, n);

    acb_init(a);
    acb_log1p(a, f, prec);

    if (flen == 1)
    {
        _acb_vec_zero(res + 1, n - 1);
    }
    else if (n == 2)
    {
        acb_add_ui(res, f + 0, 1, prec);
        acb_div(res + 1, f + 1, res + 0, prec);
    }
    else if (_acb_vec_is_zero(f + 1, flen - 2))  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        acb_add_ui(res, f + 0, 1, prec);

        for (i = 1, j = d; j < n; j += d, i++)
        {
            if (i == 1)
                acb_div(res + j, f + d, res, prec);
            else
                acb_mul(res + j, res + j - d, res + d, prec);
            _acb_vec_zero(res + j - d + 1, flen - 2);
        }
        _acb_vec_zero(res + j - d + 1, n - (j - d + 1));

        for (i = 2, j = 2 * d; j < n; j += d, i++)
            acb_div_si(res + j, res + j, i % 2 ? i : -i, prec);
    }
    else
    {
        acb_ptr f_diff, f_inv;
        slong alloc;

        alloc = n + flen;
        f_inv = _acb_vec_init(alloc);
        f_diff = f_inv + n;

        acb_add_ui(f_diff, f, 1, prec);
        _acb_vec_set(f_diff + 1, f + 1, flen - 1);
        _acb_poly_inv_series(f_inv, f_diff, flen, n, prec);
        _acb_poly_derivative(f_diff, f, flen, prec);
        _acb_poly_mullow(res, f_inv, n - 1, f_diff, flen - 1, n - 1, prec);
        _acb_poly_integral(res, res, n, prec);

        _acb_vec_clear(f_inv, alloc);
    }

    acb_swap(res, a);
    acb_clear(a);
}

void
acb_poly_log1p_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    slong flen = f->length;

    if (flen == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    if (flen == 1 /*&& !acb_contains_si(f->coeffs, -1)*/)
        n = 1;

    acb_poly_fit_length(res, n);
    _acb_poly_log1p_series(res->coeffs, f->coeffs, flen, n, prec);
    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

