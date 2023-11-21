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
_acb_hypgeom_shi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec)
{
    acb_t c;

    acb_init(c);
    acb_hypgeom_shi(c, h, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_mul_onei(g, h);
        acb_sinc(g, g, prec);
        acb_mul(g + 1, g, h + 1, prec);
    }
    else
    {
        acb_ptr t, u;
        slong i;

        t = _acb_vec_init(len - 1);
        u = _acb_vec_init(hlen);

        /* Shi(h(x)) = integral(h'(x) sinc(i h(x))) */
        for (i = 0; i < hlen; i++)
            acb_mul_onei(u + i, h + i);
        _acb_poly_sinc_series(t, u, hlen, len - 1, prec);
        _acb_poly_derivative(u, h, hlen, prec);
        _acb_poly_mullow(g, t, len - 1, u, hlen - 1, len - 1, prec);
        _acb_poly_integral(g, g, len, prec);

        _acb_vec_clear(t, len - 1);
        _acb_vec_clear(u, hlen);
    }

    acb_swap(g, c);
    acb_clear(c);
}

void
acb_hypgeom_shi_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        acb_poly_zero(g);
        return;
    }

    acb_poly_fit_length(g, len);
    _acb_hypgeom_shi_series(g->coeffs, h->coeffs, hlen, len, prec);
    _acb_poly_set_length(g, len);
    _acb_poly_normalise(g);
}

