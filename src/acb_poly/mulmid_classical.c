/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_mulmid_classical(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong nlo, slong nhi, slong prec)
{
    slong len = nhi - nlo;

    if (nhi == 1)
    {
        acb_mul(res, poly1, poly2, prec);
    }
    else if (poly1 == poly2 && len1 == len2)
    {
        slong i, start, stop;

        if (nlo == 0)
            acb_sqr(res, poly1, prec);

        for (i = FLINT_MAX(nlo, 1); i < FLINT_MIN(nhi, 2 * len1 - 2); i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            acb_dot(res + i - nlo, NULL, 0, poly1 + start, 1,
                poly1 + i - start, -1, stop - start + 1, prec);
            acb_mul_2exp_si(res + i - nlo, res + i - nlo, 1);
            if (i % 2 == 0 && i / 2 < len1)
                acb_addmul(res + i - nlo, poly1 + i / 2, poly1 + i / 2, prec);
        }

        if (nhi >= 2 * len1 - 1)
            acb_sqr(res + 2 * len1 - 2 - nlo, poly1 + len1 - 1, prec);
    }
    else if (len1 == 1)
    {
        _acb_vec_scalar_mul(res, poly2 + nlo, len, poly1, prec);
    }
    else if (len2 == 1)
    {
        _acb_vec_scalar_mul(res, poly1 + nlo, len, poly2, prec);
    }
    else
    {
        slong i, top1, top2;

        for (i = nlo; i < nhi; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            acb_dot(res + i - nlo, NULL, 0, poly1 + i - top2, 1,
                poly2 + top2, -1, top1 + top2 - i + 1, prec);
        }
    }
}

void
acb_poly_mulmid_classical(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong nlo, slong nhi, slong prec)
{
    slong xlen, ylen, zlen;

    xlen = poly1->length;
    ylen = poly2->length;

    if (xlen == 0 || ylen == 0 || nlo >= FLINT_MIN(nhi, xlen + ylen - 1))
    {
        acb_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, xlen + ylen - 1);
    zlen = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        acb_poly_t tmp;
        acb_poly_init2(tmp, zlen);
        _acb_poly_mulmid_classical(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
        acb_poly_swap(res, tmp);
        acb_poly_clear(tmp);
    }
    else
    {
        acb_poly_fit_length(res, zlen);
        _acb_poly_mulmid_classical(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
    }

    _acb_poly_set_length(res, zlen);
    _acb_poly_normalise(res);
}
