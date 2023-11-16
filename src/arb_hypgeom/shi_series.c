/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

/* todo: use a sinch function? */
void
_arb_hypgeom_shi_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        arb_hypgeom_shi(g, h, prec);
        _arb_vec_zero(g + 1, len - 1);
    }
    else
    {
        acb_ptr t;
        slong i;

        t = _acb_vec_init(len);

        for (i = 0; i < hlen; i++)
            arb_set(acb_realref(t + i), h + i);

        _acb_hypgeom_shi_series(t, t, hlen, len, prec);

        for (i = 0; i < len; i++)
            arb_swap(g + i, acb_realref(t + i));

        _acb_vec_clear(t, len);
    }
}

void
arb_hypgeom_shi_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, len);
    _arb_hypgeom_shi_series(g->coeffs, h->coeffs, hlen, len, prec);
    _arb_poly_set_length(g, len);
    _arb_poly_normalise(g);
}

