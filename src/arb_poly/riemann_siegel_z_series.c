/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
_arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec)
{
    slong i, alloc;
    arb_ptr t, u, v, w, q;

    hlen = FLINT_MIN(hlen, len);

    alloc = 5 * len;
    t = _arb_vec_init(alloc);
    u = t + len;
    v = u + len;
    w = v + len;
    q = w + len;

    /* (v + wi) = zeta(1/2 + i (s+x)) */
    {
        acb_ptr sx, z;
        acb_t a;
        slong slen = FLINT_MIN(len, 2);

        z = _acb_vec_init(len);
        sx = _acb_vec_init(slen);
        acb_init(a);

        acb_one(a);
        acb_one(sx);
        acb_mul_2exp_si(sx, sx, -1);
        arb_set(acb_imagref(sx), h);
        if (slen > 1)
            arb_one(acb_imagref(sx + 1));

        _acb_poly_zeta_series(z, sx, slen, a, 0, len, prec);

        for (i = 0; i < len; i++)
        {
            arb_set(v + i, acb_realref(z + i));
            arb_set(w + i, acb_imagref(z + i));
        }

        acb_clear(a);
        _acb_vec_clear(z, len);
        _acb_vec_clear(sx, slen);
    }

    /* (t + ui) = exp(i theta(s+x)) */
    arb_set(u, h);
    if (len > 1)
        arb_one(u + 1);
    _arb_poly_riemann_siegel_theta_series(t, u, 2, len, prec);
    _arb_poly_sin_cos_series(u, t, t, len, len, prec);

    /* (t+ui)(v+wi) = (tv-uw) + (tw+uv)i */
    _arb_poly_mullow(q, t, len, v, len, len, prec);
    _arb_poly_mullow(t, u, len, w, len, len, prec);
    _arb_vec_sub(t, q, t, len, prec);

    /* compose with nonconstant part */
    arb_zero(u);
    _arb_vec_set(u + 1, h + 1, hlen - 1);
    _arb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _arb_vec_clear(t, alloc);
}

void
arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_fit_length(res, n);

    if (f->length == 0)
    {
        arb_t t;
        arb_init(t);
        _arb_poly_riemann_siegel_z_series(res->coeffs, t, 1, n, prec);
        arb_clear(t);
    }
    else
    {
        _arb_poly_riemann_siegel_z_series(res->coeffs, f->coeffs, f->length, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

