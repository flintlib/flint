/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re)
{
    slong i, len;

    len = arb_poly_length(re);

    acb_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        arb_set(acb_realref(poly->coeffs + i), re->coeffs + i);
        arb_zero(acb_imagref(poly->coeffs + i));
    }

    _acb_poly_set_length(poly, len);
}

void
acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im)
{
    slong i, rlen, ilen, len;

    rlen = arb_poly_length(re);
    ilen = arb_poly_length(im);
    len = FLINT_MAX(rlen, ilen);

    acb_poly_fit_length(poly, len);

    for (i = 0; i < rlen; i++)
        arb_set(acb_realref(poly->coeffs + i), re->coeffs + i);
    for (i = rlen; i < len; i++)
        arb_zero(acb_realref(poly->coeffs + i));

    for (i = 0; i < ilen; i++)
        arb_set(acb_imagref(poly->coeffs + i), im->coeffs + i);
    for (i = ilen; i < len; i++)
        arb_zero(acb_imagref(poly->coeffs + i));

    _acb_poly_set_length(poly, len);
}
