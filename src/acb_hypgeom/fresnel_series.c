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
_acb_hypgeom_fresnel_series(acb_ptr s, acb_ptr c,
    acb_srcptr h, slong hlen, int normalized, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_hypgeom_fresnel(s, c, h, normalized, prec);

        if (s != NULL) _acb_vec_zero(s + 1, len - 1);
        if (c != NULL) _acb_vec_zero(c + 1, len - 1);
    }
    else
    {
        acb_t s0, c0;
        acb_ptr t, u, v;
        slong ulen;

        acb_init(s0);
        acb_init(c0);

        acb_hypgeom_fresnel((s != NULL) ? s0 : NULL,
                            (c != NULL) ? c0 : NULL, h, normalized, prec);

        t = _acb_vec_init(len);
        u = _acb_vec_init(len);
        v = _acb_vec_init(len);

        /* normalized: */
        /* C(h(x)) = integral(h'(x) cos(-(pi/2) h(x)^2)) */
        /* S(h(x)) = -integral(h'(x) sin(-(pi/2) h(x)^2)) */
        ulen = FLINT_MIN(len, 2 * hlen - 1);

        _acb_poly_mullow(u, h, hlen, h, hlen, ulen, prec);
        _acb_vec_neg(u, u, ulen);

        if (normalized)
        {
            _acb_vec_scalar_mul_2exp_si(u, u, ulen, -1);
            _acb_poly_sin_cos_pi_series(u, v, u, ulen, len, prec);
        }
        else
        {
            _acb_poly_sin_cos_series(u, v, u, ulen, len, prec);
        }

        _acb_poly_derivative(t, h, hlen, prec);

        if (s != NULL)
        {
            _acb_poly_mullow(s, u, len, t, hlen - 1, len, prec);
            _acb_poly_integral(s, s, len, prec);
            _acb_vec_neg(s, s, len);
            acb_swap(s, s0);
        }

        if (c != NULL)
        {
            _acb_poly_mullow(c, v, len, t, hlen - 1, len, prec);
            _acb_poly_integral(c, c, len, prec);
            acb_swap(c, c0);
        }

        _acb_vec_clear(t, len);
        _acb_vec_clear(u, len);
        _acb_vec_clear(v, len);

        acb_clear(s0);
        acb_clear(c0);
    }
}

void
acb_hypgeom_fresnel_series(acb_poly_t s, acb_poly_t c,
    const acb_poly_t h, int normalized, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        if (s != NULL) acb_poly_zero(s);
        if (c != NULL) acb_poly_zero(c);
        return;
    }

    if (s != NULL) acb_poly_fit_length(s, len);
    if (c != NULL) acb_poly_fit_length(c, len);

    _acb_hypgeom_fresnel_series((s != NULL) ? s->coeffs : NULL,
                                (c != NULL) ? c->coeffs : NULL,
                        h->coeffs, hlen, normalized, len, prec);

    if (s != NULL) _acb_poly_set_length(s, len);
    if (c != NULL) _acb_poly_set_length(c, len);
    if (s != NULL) _acb_poly_normalise(s);
    if (c != NULL) _acb_poly_normalise(c);
}

