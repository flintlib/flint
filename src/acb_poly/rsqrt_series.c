/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
_acb_poly_rsqrt_series(acb_ptr g,
    acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    while (hlen > 0 && acb_is_zero(h + hlen - 1))
        hlen--;

    if (hlen <= 1)
    {
        acb_rsqrt(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_rsqrt(g, h, prec);
        acb_div(g + 1, h + 1, h, prec);
        acb_mul(g + 1, g + 1, g, prec);
        acb_mul_2exp_si(g + 1, g + 1, -1);
        acb_neg(g + 1, g + 1);
    }
    else if (_acb_vec_is_zero(h + 1, hlen - 2))
    {
        acb_t t;
        acb_init(t);
        arf_set_si_2exp_si(arb_midref(acb_realref(t)), -1, -1);
        _acb_poly_binomial_pow_acb_series(g, h, hlen, t, len, prec);
        acb_clear(t);
    }
    else
    {
        /* Third-order Newton iteration: with e = 1 - h g^2 = O(x^m),
           h^(-1/2) = g (1 + e/2 + 3e^2/8) + O(x^(3m)). */
        acb_ptr t, u, v;
        slong a[FLINT_BITS];
        slong i, m, n, k, tlen, hnlen, plen;

        t = _acb_vec_init(2 * len + (len + 2) / 3);
        u = t + len;
        v = u + len;

        acb_rsqrt(g, h, prec);

        a[i = 0] = n = len;
        while (n > 1)
            a[++i] = (n = (n + 2) / 3);

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];
            k = n - 2 * m;

            tlen = FLINT_MIN(2 * m - 1, n);
            _acb_poly_mullow(t, g, m, g, m, tlen, prec);

            /* u = (h g^2)[m, n) = -e / x^m */
            hnlen = FLINT_MIN(hlen, n);
            plen = FLINT_MIN(n, tlen + hnlen - 1);
            _acb_poly_mulmid(u, t, tlen, h, hnlen, m, plen, prec);
            _acb_vec_zero(u + plen - m, n - plen);

            if (k > 0)
            {
                /* e/2 + 3e^2/8 = -(4u - 3 x^m u^2) / 8 */
                _acb_poly_mullow(v, u, k, u, k, k, prec);
                _acb_vec_scalar_mul_si(v, v, k, 3, prec);
                _acb_vec_scalar_mul_2exp_si(u, u, n - m, 2);
                _acb_vec_sub(u + m, u + m, v, k, prec);
            }

            _acb_poly_mullow(g + m, g, FLINT_MIN(m, n - m), u, n - m, n - m, prec);
            _acb_vec_scalar_mul_2exp_si(g + m, g + m, n - m, (k > 0) ? -3 : -1);
            _acb_vec_neg(g + m, g + m, n - m);
        }

        _acb_vec_clear(t, 2 * len + (len + 2) / 3);
    }
}

void
acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        acb_poly_t t;
        acb_poly_init(t);
        acb_poly_rsqrt_series(t, h, n, prec);
        acb_poly_swap(g, t);
        acb_poly_clear(t);
        return;
    }

    acb_poly_fit_length(g, n);
    if (h->length == 0)
        _acb_vec_indeterminate(g->coeffs, n);
    else
        _acb_poly_rsqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
