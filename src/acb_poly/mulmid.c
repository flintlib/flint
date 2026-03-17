/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb_poly.h"

void
_acb_poly_mulmid(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong nlo, slong nhi, slong prec)
{
    if (len1 <= 7 || len2 <= 7 || nhi <= 7)
    {
        _acb_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, prec);
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

        /* todo: tuning copied from mullow; needs retuning for small nhi - nlo */
        if (2 * FLINT_MIN(len1, len2) <= cutoff || nhi <= cutoff)
            _acb_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, prec);
        else
            _acb_poly_mulmid_transpose(res, poly1, len1, poly2, len2, nlo, nhi, prec);
    }
}

void
acb_poly_mulmid(acb_poly_t res, const acb_poly_t poly1,
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
        _acb_poly_mulmid(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
        acb_poly_swap(res, tmp);
        acb_poly_clear(tmp);
    }
    else
    {
        acb_poly_fit_length(res, zlen);
        _acb_poly_mulmid(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
    }

    _acb_poly_set_length(res, zlen);
    _acb_poly_normalise(res);
}
