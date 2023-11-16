/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

void
_arb_hypgeom_erfc_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
{
    arb_t c;
    arb_init(c);

    arb_hypgeom_erfc(c, h, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, len - 1);
    }
    else
    {
        arb_ptr t, u;
        slong ulen;

        t = _arb_vec_init(len);
        u = _arb_vec_init(len);

        /* erfc(h(x)) = integral(-h'(x) exp(-h(x)^2)) * (2/sqrt(pi)) */
        ulen = FLINT_MIN(len, 2 * hlen - 1);

        _arb_poly_mullow(u, h, hlen, h, hlen, ulen, prec);
        _arb_vec_neg(u, u, ulen);
        _arb_poly_exp_series(u, u, ulen, len, prec);
        _arb_poly_derivative(t, h, hlen, prec);
        _arb_poly_mullow(g, u, len, t, hlen - 1, len, prec);
        _arb_poly_integral(g, g, len, prec);

        arb_const_sqrt_pi(t, prec);
        arb_inv(t, t, prec);
        arb_mul_2exp_si(t, t, 1);
        _arb_vec_scalar_mul(g, g, len, t, prec);
        _arb_vec_neg(g, g, len);

        _arb_vec_clear(t, len);
        _arb_vec_clear(u, len);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_hypgeom_erfc_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (len == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(g);
        return;
    }

    arb_poly_fit_length(g, len);
    _arb_hypgeom_erfc_series(g->coeffs, h->coeffs, hlen, len, prec);
    _arb_poly_set_length(g, len);
    _arb_poly_normalise(g);
}

