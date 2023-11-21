/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_shift_left(acb_ptr res, acb_srcptr poly, slong len, slong n)
{
    slong i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = len; i--; )
            acb_set(res + n + i, poly + i);
    }
    else
    {
        for (i = len; i--; )
            acb_swap(res + n + i, res + i);
    }

    for (i = 0; i < n; i++)
        acb_zero(res + i);
}

void
acb_poly_shift_left(acb_poly_t res, const acb_poly_t poly, slong n)
{
    if (n == 0)
    {
        acb_poly_set(res, poly);
        return;
    }

    if (poly->length == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, poly->length + n);
    _acb_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n);
    _acb_poly_set_length(res, poly->length + n);
}

