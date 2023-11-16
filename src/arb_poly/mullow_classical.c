/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_mullow_classical(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec)
{
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        arb_mul(res, poly1, poly2, prec);
    }
    else if (poly1 == poly2 && len1 == len2)
    {
        slong i, start, stop;

        arb_sqr(res, poly1, prec);
        arb_mul(res + 1, poly1, poly1 + 1, prec);
        arb_mul_2exp_si(res + 1, res + 1, 1);

        for (i = 2; i < FLINT_MIN(n, 2 * len1 - 3); i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            arb_dot(res + i, NULL, 0, poly1 + start, 1,
                poly1 + i - start, -1, stop - start + 1, prec);
            arb_mul_2exp_si(res + i, res + i, 1);
            if (i % 2 == 0 && i / 2 < len1)
                arb_addmul(res + i, poly1 + i / 2, poly1 + i / 2, prec);
        }

        if (len1 > 2 && n >= 2 * len1 - 2)
        {
            arb_mul(res + 2 * len1 - 3, poly1 + len1 - 1, poly1 + len1 - 2, prec);
            arb_mul_2exp_si(res + 2 * len1 - 3, res + 2 * len1 - 3, 1);
        }

        if (n >= 2 * len1 - 1)
            arb_sqr(res + 2 * len1 - 2, poly1 + len1 - 1, prec);
    }
    else if (len1 == 1)
    {
        _arb_vec_scalar_mul(res, poly2, n, poly1, prec);
    }
    else if (len2 == 1)
    {
        _arb_vec_scalar_mul(res, poly1, n, poly2, prec);
    }
    else
    {
        slong i, top1, top2;

        arb_mul(res, poly1, poly2, prec);

        for (i = 1; i < n; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            arb_dot(res + i, NULL, 0, poly1 + i - top2, 1,
                poly2 + top2, -1, top1 + top2 - i + 1, prec);
        }
    }
}

void
arb_poly_mullow_classical(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                slong n, slong prec)
{
    slong len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_mullow_classical(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_mullow_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}
