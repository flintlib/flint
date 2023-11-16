/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

void
_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (n == 1)
    {
        acb_mul(res, poly1, poly2, prec);
    }
    else if (n <= 7 || len1 <= 7 || len2 <= 7)
    {
        _acb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
    }
    else
    {
        slong cutoff;
        double p;

        if (prec <= 2 * FLINT_BITS)
        {
            cutoff = 110;
        }
        else
        {
            p = log(prec);

            cutoff = 10000.0 / (p * p * p);
            cutoff = FLINT_MIN(cutoff, 60);
            if (poly1 == poly2 && prec >= 256)
                cutoff *= 1.25;
            if (poly1 == poly2 && prec >= 4096)
                cutoff *= 1.25;
            cutoff = FLINT_MAX(cutoff, 8);
        }

        if (2 * FLINT_MIN(len1, len2) <= cutoff || n <= cutoff)
            _acb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else
            _acb_poly_mullow_transpose(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
acb_poly_mullow(acb_poly_t res, const acb_poly_t poly1,
                    const acb_poly_t poly2, slong n, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    /* Hack to avoid temporary allocations with first derivatives. */
    if (n <= 2 && !(len1 == 2 && len2 == 2))
    {
        acb_poly_fit_length(res, n);

        if (n == 1)
        {
            acb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len2 == 1)
        {
            acb_mul(res->coeffs + 1, poly1->coeffs + 1, poly2->coeffs, prec);
            acb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len1 == 1)
        {
            acb_mul(res->coeffs + 1, poly2->coeffs + 1, poly1->coeffs, prec);
            acb_mul(res->coeffs, poly2->coeffs, poly1->coeffs, prec);
        }
        else
        {
            if (res == poly1 || res == poly2)
            {
                acb_t t;
                acb_init(t);
                acb_mul(t, poly1->coeffs, poly2->coeffs + 1, prec);
                acb_addmul(t, poly2->coeffs, poly1->coeffs + 1, prec);
                acb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                acb_swap(t, res->coeffs + 1);
                acb_clear(t);
            }
            else
            {
                acb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                acb_mul(res->coeffs + 1, poly1->coeffs, poly2->coeffs + 1, prec);
                acb_addmul(res->coeffs + 1, poly2->coeffs, poly1->coeffs + 1, prec);
            }
        }

        _acb_poly_set_length(res, n);
        _acb_poly_normalise(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_mullow(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_mullow(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

