/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

/* Compute R_C(1,1+x) = atan(sqrt(x))/sqrt(x) = 2F1(1,1/2,3/2,-x) */
static void
_acb_elliptic_rc1(acb_t res, const acb_t x, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_sqrt(t, x, prec + 2);
    acb_atan(res, t, prec + 2);
    acb_div(res, res, t, prec);
    acb_clear(t);
}

void
acb_elliptic_rc1(acb_t res, const acb_t x, slong prec)
{
    mag_t xm;

    mag_init(xm);
    acb_get_mag(xm, x);

    if (mag_cmp_2exp_si(xm, 0) < 0)
    {
        slong k, n;

        for (n = 1; n <= 6; n++)
        {
            if (mag_cmp_2exp_si(xm, -prec / n) < 0)
                break;
        }

        /* Use Taylor series: 1 - x/3 + x^2/5 - x^3/7 + x^4/9 + ... */
        if (n <= 6)
        {
            const short coeffs[] = {3465, -1155, 693, -495, 385, -315};

            acb_t s;
            acb_init(s);

            for (k = n - 1; k >= 0; k--)
            {
                acb_mul(s, s, x, prec);
                acb_add_si(s, s, coeffs[k], prec);
            }

            acb_div_si(s, s, coeffs[0], prec);
            mag_geom_series(xm, xm, n);
            if (acb_is_real(x))
                arb_add_error_mag(acb_realref(s), xm);
            else
                acb_add_error_mag(s, xm);

            acb_set(res, s);
            acb_clear(s);
        }
        else if (acb_is_exact(x))
        {
            _acb_elliptic_rc1(res, x, prec);
        }
        else
        {
            acb_t w;
            mag_t err, rad;

            acb_init(w);
            mag_init(err);
            mag_init(rad);

            /* On the unit disc, |f'(x)| <= 0.5 / |1+x| */
            acb_add_ui(w, x, 1, MAG_BITS);
            acb_get_mag_lower(err, w);
            mag_one(rad);
            mag_mul_2exp_si(rad, rad, -1);
            mag_div(err, rad, err);
            mag_hypot(rad, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
            mag_mul(err, err, rad);

            acb_set(w, x);
            mag_zero(arb_radref(acb_realref(w)));
            mag_zero(arb_radref(acb_imagref(w)));
            _acb_elliptic_rc1(w, w, prec);

            if (acb_is_real(x))
                arb_add_error_mag(acb_realref(w), err);
            else
                acb_add_error_mag(w, err);

            acb_set(res, w);

            acb_clear(w);
            mag_clear(err);
            mag_clear(rad);
        }
    }
    else
    {
        _acb_elliptic_rc1(res, x, prec);
    }

    mag_clear(xm);
}

