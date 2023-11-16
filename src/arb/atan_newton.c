/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_fmpz_poly.h"

static const fmpz atan_coeffs[] =
{
    334639305, -111546435, 66927861, -47805615, 37182145,
    -30421755, 25741485, -22309287, 19684665, -17612595,
    15935205, -14549535
};

static const ulong atan_den = 334639305;

void
_arb_atan_taylor(arb_t res, const arb_t x, slong prec)
{
    slong mag;
    slong n;
    mag_t err;
    arb_t t;

    mag = arf_abs_bound_lt_2exp_si(arb_midref(x));

    if (mag > -1)
    {
        arb_indeterminate(res);
        return;
    }

    arb_init(t);
    mag_init(err);

    if (mag < -prec)
        n = 1;
    else
        n = (prec + 2 * (-mag) - 1) / (2 * (-mag));

    n = FLINT_MIN(n, 12);

    arb_get_mag(err, x);
    mag_geom_series(err, err, 2 * n + 1);
    mag_div_ui(err, err, 2 * n + 1);

    arb_mul(t, x, x, prec + 10);
    _arb_fmpz_poly_evaluate_arb_rectangular(t, atan_coeffs, n, t, prec + 10);
    arb_mul(res, t, x, prec + 10);
    arb_div_ui(res, res, atan_den, prec);

    arb_add_error_mag(res, err);

    mag_clear(err);
    arb_clear(t);
}

void
arb_atan_newton(arb_t res, const arb_t x, slong prec)
{
    arb_t t, s, c, w;
    slong n, wp, wp2;
    mag_t err;
    slong xmag, extra;

    if (arb_is_zero(x))
    {
        arb_zero(res);
        return;
    }

    if (!arb_is_finite(x))
    {
        arb_indeterminate(res);
        return;
    }

    xmag = arf_abs_bound_lt_2exp_si(arb_midref(x));

    if (xmag > 4)
    {
        int sgn = arf_sgn(arb_midref(x));

        if (arb_contains_zero(x))
        {
            arb_indeterminate(res);
            return;
        }

        wp = FLINT_MAX(0, prec - xmag) + 15;

        arb_init(t);
        arb_inv(t, x, wp);
        arb_atan_newton(t, t, wp);
        arb_const_pi(res, prec + 15);
        arb_mul_2exp_si(res, res, -1);
        if (sgn < 0)
            arb_neg(res, res);
        arb_sub(res, res, t, prec);
        arb_clear(t);
        return;
    }

    arb_init(t);
    arb_init(s);
    arb_init(c);
    arb_init(w);
    mag_init(err);

    /* quick and accurate calculation near 0 */
    if (xmag < -prec / 20)
    {
        _arb_atan_taylor(res, x, prec);
    }
    /* prec < ARB_ATAN_NEWTON_PREC / 2 if we recursed to this function,
       but we actually recurse by calling arb_atan, so this just
       needs to be a failsafe */
    else if (prec <= 64)
    {
        arb_atan(res, x, prec);
    }
    else
    {
        if (prec <= 6000)
            n = 5;
        else if (prec <= 100000)
            n = 7;
        else if (prec <= 1000000)
            n = 9;
        else
            n = 11;

        wp = prec + 10 + (-xmag);
        wp2 = wp * (n - 1) / n;

        extra = 10;

        arb_atan(t, x, wp / n + extra);
        mag_zero(arb_radref(t));

        arb_sin_cos(s, c, t, wp);
        arb_set_round(res, x, wp);

        arb_mul(w, c, res, wp);
        arb_sub(w, w, s, wp2);
        arb_mul(res, s, res, wp);
        arb_add(res, res, c, wp2);
        arb_div(w, w, res, wp2);

        _arb_atan_taylor(res, w, wp2);

        arb_add(res, t, res, prec);
    }

    arb_clear(t);
    arb_clear(s);
    arb_clear(c);
    arb_clear(w);
    mag_clear(err);
}

void
arb_atan_arf_newton(arb_t res, const arf_t x, slong prec)
{
    arb_set_arf(res, x);
    arb_atan_newton(res, res, prec);
}
