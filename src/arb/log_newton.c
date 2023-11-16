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

static const fmpz log_coeffs[] =
{
    0, 232792560, -116396280, 77597520, -58198140, 46558512,
    -38798760, 33256080, -29099070, 25865840, -23279256, 21162960,
    -19399380, 17907120, -16628040, 15519504, -14549535, 13693680,
    -12932920, 12252240, -11639628
};

static const ulong log_den = 232792560;

void
arb_log_newton(arb_t res, const arb_t x, slong prec)
{
    arb_t t, w;
    slong wp, wp2;
    mag_t err;
    slong n, x1mag, extra, ebits;

    if (arb_is_one(x))
    {
        arb_zero(res);
        return;
    }

    arb_init(t);
    arb_init(w);
    mag_init(err);

    arf_sub_ui(arb_midref(t), arb_midref(x), 1, 8, ARF_RND_DOWN);
    x1mag = arf_abs_bound_lt_2exp_si(arb_midref(t));

    /* quick and accurate calculation near 1 */
    if (x1mag < -prec / 16)
    {
        n = prec / (-x1mag) + 2;

        arb_sub_ui(t, x, 1, prec + 10);

        arb_get_mag(err, t);
        mag_geom_series(err, err, n);
        mag_div_ui(err, err, n);

        _arb_fmpz_poly_evaluate_arb_rectangular(res, log_coeffs, n, t, prec + 10);
        arb_div_ui(res, res, log_den, prec);

        arb_add_error_mag(res, err);
    }
    /* catch well below ARB_LOG_NEWTON_PREC so that we never hit
       infinite recursion */
    else if (prec <= ARB_LOG_NEWTON_PREC / 2)
    {
        arb_log(res, x, prec);
    }
    else
    {
        if (prec <= 3200)
            n = 4;
        else if (prec <= 6000)
            n = 6;
        else if (prec <= 300000)
            n = 7;
        else if (prec <= 1000000)
            n = 9;
        else
            n = 12;

        n = FLINT_MAX(n, 2);
        n = FLINT_MIN(n, 16);

        extra = 10;
        ebits = fmpz_bits(ARF_EXPREF(arb_midref(x)));
        extra += ebits;

        if (extra > 30)
        {
            fmpz_t q;
            fmpz_init(q);
            fmpz_set(q, ARF_EXPREF(arb_midref(x)));
            fmpz_neg(q, q);
            arb_mul_2exp_fmpz(t, x, q);
            arb_log_newton(res, t, prec + 5 - ebits);
            arb_const_log2(t, prec + 5);
            arb_submul_fmpz(res, t, q, prec);
            fmpz_clear(q);
        }
        else
        {
            wp = prec + 10;
            if (x1mag < 0)
                wp += (-x1mag);

            wp2 = wp * (n - 1) / n;

            arb_log(t, x, wp / n + extra);
            mag_zero(arb_radref(t));

            arb_neg(w, t);
            arb_exp(w, w, wp);
            arb_set_round(res, x, wp);
            arb_mul(w, w, res, wp);
            arb_sub_ui(w, w, 1, wp2);

            arb_get_mag(err, w);
            mag_geom_series(err, err, n);
            mag_div_ui(err, err, n);

            _arb_fmpz_poly_evaluate_arb_rectangular(res, log_coeffs, n, w, wp2);
            arb_div_ui(res, res, log_den, wp2);
            arb_add_error_mag(res, err);

            arb_add(res, t, res, prec);
        }
    }

    arb_clear(t);
    arb_clear(w);
    mag_clear(err);
}

void
arb_log_arf_newton(arb_t res, const arf_t x, slong prec)
{
    arb_set_arf(res, x);
    arb_log_newton(res, res, prec);
}
