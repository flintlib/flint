/*
    Copyright (C) 2016, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "acb_poly.h"

void
_acb_sinc_jet_zero(acb_ptr res, const acb_t z, slong len, slong prec)
{
    mag_ptr D;
    slong n;
    int is_real, is_imag;
    mag_t zm, err, fac;
    arb_t wide;

    is_real = acb_is_real(z);
    is_imag = arb_is_zero(acb_realref(z));

    mag_init(zm);
    mag_init(err);
    mag_init(fac);
    arb_init(wide);

    acb_get_mag(zm, z);

    D = _mag_vec_init(len + 1);

    /* Compute sequence of derivative bounds

        sinc^{(n)}(z) = int_0^1 t^n cos(t z + n pi / 2) dt

        |sinc^{(n)}(z)| <= cosh(im(z)) min(1/|z|, 1/(n+1))
        |sinc^{(n)}(z)| <= cosh(im(z)) min(1/|z|, 1/(n+1), |z|/(n+2))  (odd n)
    */
    {
        mag_t C, b1, b2, b3;

        mag_init(C);
        mag_init(b1);
        mag_init(b2);
        mag_init(b3);

        /* C = cosh(im(z)) */
        arb_get_mag(C, acb_imagref(z));
        mag_cosh(C, C);

        /* b1 = 1/|z| */
        acb_get_mag_lower(b1, z);
        mag_inv(b1, b1);

        for (n = 0; n <= len; n++)
        {
            /* b2 = 1/(n+1) */
            mag_one(b2);
            mag_div_ui(b2, b2, n + 1);
            mag_min(b2, b1, b2);

            /* b3 = |z|/(n+2) */
            if (n % 2 == 1)
            {
                mag_div_ui(b3, zm, n + 2);
                mag_min(b2, b2, b3);
            }

            mag_mul(D + n, C, b2);
        }

        mag_clear(C);
        mag_clear(b1);
        mag_clear(b2);
        mag_clear(b3);
    }

    mag_one(fac);

    for (n = 0; n < len; n++)
    {
        /* sinc^{(n)}(0) / n! */
        if (n == 0)
            acb_one(res + n);
        else if (n % 2 == 1)
            acb_zero(res + n);
        else
            acb_div_si(res + n, res + n - 2, -n * (n + 1), prec);

        if (n > 1)
            mag_div_ui(fac, fac, n);

        /* |sinc^{(n)}(0 + eps) - sinc^{(n)}(0)| / n! <= |eps| * |sinc^{(n+1)}(z)| / n!, z = (+/- eps) */
        mag_mul(err, zm, D + n + 1);
        mag_mul(err, err, fac);
        acb_add_error_mag(res + n, err);

        /* sinc^{(n)}(+/- eps) is a possibly better enclosure */
        arb_zero(wide);
        mag_mul(err, D + n, fac);
        arb_add_error_mag(wide, err);
        arb_intersection(acb_realref(res + n), acb_realref(res + n), wide, prec);
        arb_intersection(acb_imagref(res + n), acb_imagref(res + n), wide, prec);

        if (is_real || (is_imag && (n % 2 == 0)))
            arb_zero(acb_imagref(res + n));

        if (is_imag && (n % 2 == 1))
            arb_zero(acb_realref(res + n));
    }

    _mag_vec_clear(D, len + 1);

    mag_clear(zm);
    mag_clear(err);
    mag_clear(fac);
    arb_clear(wide);
}

void
_acb_poly_sinc_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sinc(g, h, prec);
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        acb_ptr t, u;

        t = _acb_vec_init(n + 1);
        u = _acb_vec_init(hlen);

        _acb_vec_set(u, h, hlen);

        if (acb_is_zero(h))
        {
            _acb_poly_sin_series(t, u, hlen, n + 1, prec);
            _acb_poly_div_series(g, t + 1, n, u + 1, hlen - 1, n, prec);
        }
        else if (acb_contains_zero(h))
        {
            _acb_sinc_jet_zero(t, h, n, prec);
            /* compose with nonconstant part */
            acb_zero(u);
            _acb_poly_compose_series(g, t, n, u, hlen, n, prec);
        }
        else
        {
            _acb_poly_sin_series(t, u, hlen, n, prec);
            _acb_poly_div_series(g, t, n, u, hlen, n, prec);
        }

        _acb_vec_clear(t, n + 1);
        _acb_vec_clear(u, hlen);
    }
}

void
acb_poly_sinc_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(g, n);
    _acb_poly_sinc_series(g->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
