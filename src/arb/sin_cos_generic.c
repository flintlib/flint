/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
_arb_sin_cos_generic(arb_t s, arb_t c, const arf_t x, const mag_t xrad, slong prec)
{
    int want_sin, want_cos;
    slong maglim;

    want_sin = (s != NULL);
    want_cos = (c != NULL);

    if (arf_is_zero(x) && mag_is_zero(xrad))
    {
        if (want_sin) arb_zero(s);
        if (want_cos) arb_one(c);
        return;
    }

    if (!arf_is_finite(x) || !mag_is_finite(xrad))
    {
        if (arf_is_nan(x))
        {
            if (want_sin) arb_indeterminate(s);
            if (want_cos) arb_indeterminate(c);
        }
        else
        {
            if (want_sin) arb_zero_pm_one(s);
            if (want_cos) arb_zero_pm_one(c);
        }
        return;
    }

    maglim = FLINT_MAX(65536, 4 * prec);

    if (mag_cmp_2exp_si(xrad, -16) > 0 || arf_cmpabs_2exp_si(x, maglim) > 0)
    {
        _arb_sin_cos_wide(s, c, x, xrad, prec);
        return;
    }

    if (arf_cmpabs_2exp_si(x, -(prec/2) - 2) <= 0)
    {
        mag_t t, u, v;
        mag_init(t);
        mag_init(u);
        mag_init(v);

        arf_get_mag(t, x);
        mag_add(t, t, xrad);
        mag_mul(u, t, t);

        /* |sin(z)-z| <= z^3/6 */
        if (want_sin)
        {
            arf_set(arb_midref(s), x);
            mag_set(arb_radref(s), xrad);
            arb_set_round(s, s, prec);
            mag_mul(v, u, t);
            mag_div_ui(v, v, 6);
            arb_add_error_mag(s, v);
        }

        /* |cos(z)-1| <= z^2/2 */
        if (want_cos)
        {
            arf_one(arb_midref(c));
            mag_mul_2exp_si(arb_radref(c), u, -1);
        }

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
        return;
    }

    if (mag_is_zero(xrad))
    {
        arb_sin_cos_arf_generic(s, c, x, prec);
    }
    else
    {
        mag_t t;
        slong exp, radexp;

        mag_init_set(t, xrad);

        exp = arf_abs_bound_lt_2exp_si(x);
        radexp = MAG_EXP(xrad);
        if (radexp < MAG_MIN_LAGOM_EXP || radexp > MAG_MAX_LAGOM_EXP)
            radexp = MAG_MIN_LAGOM_EXP;

        if (want_cos && exp < -2)
            prec = FLINT_MIN(prec, 20 - FLINT_MAX(exp, radexp) - radexp);
        else
            prec = FLINT_MIN(prec, 20 - radexp);

        arb_sin_cos_arf_generic(s, c, x, prec);

        /* todo: could use quadratic bound */
        if (want_sin) mag_add(arb_radref(s), arb_radref(s), t);
        if (want_cos) mag_add(arb_radref(c), arb_radref(c), t);

        mag_clear(t);
    }
}

void
arb_sin_cos_generic(arb_t s, arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos_generic(s, c, arb_midref(x), arb_radref(x), prec);
}
