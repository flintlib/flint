/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void _acb_sinc_jet_zero(acb_ptr res, const acb_t z, slong len, slong prec);

void
_arb_sinc_jet_zero(arb_ptr res, const arb_t z, slong len, slong prec)
{
    acb_ptr t;
    slong i;
    t = _acb_vec_init(len + 1);
    acb_set_arb(t + len, z);
    _acb_sinc_jet_zero(t, t + len, len, prec);
    for (i = 0; i < len; i++)
        arb_swap(res + i, acb_realref(t + i));
    _acb_vec_clear(t, len + 1);

}

void
_arb_poly_sinc_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_sinc(g, h, prec);
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        arb_ptr t, u;

        t = _arb_vec_init(n + 1);
        u = _arb_vec_init(hlen);

        _arb_vec_set(u, h, hlen);

        if (arb_is_zero(h))
        {
            _arb_poly_sin_series(t, u, hlen, n + 1, prec);
            _arb_poly_div_series(g, t + 1, n, u + 1, hlen - 1, n, prec);
        }
        else if (arb_contains_zero(h))
        {
            _arb_sinc_jet_zero(t, h, n, prec);
            /* compose with nonconstant part */
            arb_zero(u);
            _arb_poly_compose_series(g, t, n, u, hlen, n, prec);
        }
        else
        {
            _arb_poly_sin_series(t, u, hlen, n, prec);
            _arb_poly_div_series(g, t, n, u, hlen, n, prec);
        }

        _arb_vec_clear(t, n + 1);
        _arb_vec_clear(u, hlen);
    }
}

void
arb_poly_sinc_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(g, n);
    _arb_poly_sinc_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}
