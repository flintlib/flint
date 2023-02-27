/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_log_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
{
    flen = FLINT_MIN(flen, n);

    if (flen == 1)
    {
        arb_log(res, f, prec);
        _arb_vec_zero(res + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_div(res + 1, f + 1, f + 0, prec);  /* safe since hlen >= 2 */
        arb_log(res, f, prec);
    }
    else if (_arb_vec_is_zero(f + 1, flen - 2))  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        for (i = 1, j = d; j < n; j += d, i++)
        {
            if (i == 1)
                arb_div(res + j, f + d, f + 0, prec);
            else
                arb_mul(res + j, res + j - d, res + d, prec);
            _arb_vec_zero(res + j - d + 1, flen - 2);
        }
        _arb_vec_zero(res + j - d + 1, n - (j - d + 1));

        for (i = 2, j = 2 * d; j < n; j += d, i++)
            arb_div_si(res + j, res + j, i % 2 ? i : -i, prec);

        arb_log(res, f, prec); /* done last to allow aliasing */
    }
    else
    {
        arb_ptr f_diff, f_inv;
        arb_t a;
        slong alloc;

        alloc = n + flen - 1;
        f_inv = _arb_vec_init(alloc);
        f_diff = f_inv + n;

        arb_init(a);
        arb_log(a, f, prec);

        _arb_poly_derivative(f_diff, f, flen, prec);
        _arb_poly_inv_series(f_inv, f, flen, n, prec);
        _arb_poly_mullow(res, f_inv, n - 1, f_diff, flen - 1, n - 1, prec);
        _arb_poly_integral(res, res, n, prec);
        arb_swap(res, a);

        arb_clear(a);
        _arb_vec_clear(f_inv, alloc);
    }
}

void
arb_poly_log_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_fit_length(res, n);
    if (f->length == 0)
        _arb_vec_indeterminate(res->coeffs, n);
    else
        _arb_poly_log_series(res->coeffs, f->coeffs, f->length, n, prec);
    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

