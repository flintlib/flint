/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_log_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
{
    flen = FLINT_MIN(flen, n);

    if (flen == 1)
    {
        acb_log(res, f, prec);
        _acb_vec_zero(res + 1, n - 1);
    }
    else if (n == 2)
    {
        acb_div(res + 1, f + 1, f + 0, prec);  /* safe since hlen >= 2 */
        acb_log(res, f, prec);
    }
    else if (_acb_vec_is_zero(f + 1, flen - 2))  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        for (i = 1, j = d; j < n; j += d, i++)
        {
            if (i == 1)
                acb_div(res + j, f + d, f + 0, prec);
            else
                acb_mul(res + j, res + j - d, res + d, prec);
            _acb_vec_zero(res + j - d + 1, flen - 2);
        }
        _acb_vec_zero(res + j - d + 1, n - (j - d + 1));

        for (i = 2, j = 2 * d; j < n; j += d, i++)
            acb_div_si(res + j, res + j, i % 2 ? i : -i, prec);

        acb_log(res, f, prec); /* done last to allow aliasing */
    }
    else
    {
        acb_ptr f_diff, f_inv;
        acb_t a;
        slong alloc;

        alloc = n + flen - 1;
        f_inv = _acb_vec_init(alloc);
        f_diff = f_inv + n;

        acb_init(a);
        acb_log(a, f, prec);

        _acb_poly_derivative(f_diff, f, flen, prec);
        _acb_poly_inv_series(f_inv, f, flen, n, prec);
        _acb_poly_mullow(res, f_inv, n - 1, f_diff, flen - 1, n - 1, prec);
        _acb_poly_integral(res, res, n, prec);
        acb_swap(res, a);

        acb_clear(a);
        _acb_vec_clear(f_inv, alloc);
    }
}

void
acb_poly_log_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);
    if (f->length == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_log_series(res->coeffs, f->coeffs, f->length, n, prec);
    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

