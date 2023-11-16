/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_erfi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_hypgeom_erfi(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else
    {
        slong k;
        acb_ptr t = _acb_vec_init(hlen);
        for (k = 0; k < hlen; k++)
            acb_mul_onei(t + k, h + k);
        _acb_hypgeom_erf_series(g, t, hlen, len, prec);
        for (k = 0; k < len; k++)
            acb_div_onei(g + k, g + k);
        _acb_vec_clear(t, hlen);
    }
}

void
acb_hypgeom_erfi_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        acb_poly_zero(g);
        return;
    }

    acb_poly_fit_length(g, len);
    _acb_hypgeom_erfi_series(g->coeffs, h->coeffs, hlen, len, prec);
    _acb_poly_set_length(g, len);
    _acb_poly_normalise(g);
}

