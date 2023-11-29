/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

void
_arb_poly_mullow(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (n == 1)
    {
        arb_mul(res, poly1, poly2, prec);
    }
    else if (n <= 7 || len1 <= 7 || len2 <= 7)
    {
        _arb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
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
            _arb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else
            _arb_poly_mullow_block(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1,
                        const arb_poly_t poly2, slong n, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    /* Hack to avoid temporary allocations with first derivatives. */
    if (n <= 2 && !(len1 == 2 && len2 == 2))
    {
        arb_poly_fit_length(res, n);

        if (n == 1)
        {
            arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len2 == 1)
        {
            arb_mul(res->coeffs + 1, poly1->coeffs + 1, poly2->coeffs, prec);
            arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len1 == 1)
        {
            arb_mul(res->coeffs + 1, poly2->coeffs + 1, poly1->coeffs, prec);
            arb_mul(res->coeffs, poly2->coeffs, poly1->coeffs, prec);
        }
        else
        {
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

#if 0
            if (res == poly1 || res == poly2)
            {
                arb_t t;
                arb_init(t);
                arb_mul(t, poly1->coeffs, poly2->coeffs + 1, prec);
                arb_addmul(t, poly2->coeffs, poly1->coeffs + 1, prec);
                arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                arb_swap(t, res->coeffs + 1);
                arb_clear(t);
            }
            else
            {
                arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                arb_mul(res->coeffs + 1, poly1->coeffs, poly2->coeffs + 1, prec);
                arb_addmul(res->coeffs + 1, poly2->coeffs, poly1->coeffs + 1, prec);
            }
#endif
        }

        _arb_poly_set_length(res, n);
        _arb_poly_normalise(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_mullow(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_mullow(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

