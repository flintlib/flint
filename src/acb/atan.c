/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

/* branch cuts are on (+/-i)*[1,inf] */
int
acb_atan_on_branch_cut(const acb_t z)
{
    arb_t unit;
    int result;

    if (!acb_is_finite(z))
        return 1;

    if (arb_is_nonzero(acb_realref(z)))
        return 0;

    if (arb_contains_si(acb_imagref(z), 1) || arb_contains_si(acb_imagref(z), -1))
        return 1;

    arb_init(unit);
    mag_one(arb_radref(unit));
    result = !arb_contains(unit, acb_imagref(z));
    arb_clear(unit);

    return result;
}

void
acb_atan(acb_t r, const acb_t z, slong prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_atan(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (!acb_is_finite(z))
    {
        acb_indeterminate(r);
    }
    else
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        if (acb_atan_on_branch_cut(z))
        {
            acb_mul_onei(u, z);
            acb_neg(t, u);
            acb_log1p(t, t, prec);
            acb_log1p(u, u, prec);
            acb_sub(t, t, u, prec);
            acb_mul_onei(t, t);
            acb_mul_2exp_si(r, t, -1);
        }
        else if (acb_is_exact(z))
        {
            acb_onei(t);
            acb_sub(t, t, z, prec);
            acb_div(t, z, t, prec);
            acb_mul_2exp_si(t, t, 1);
            acb_log1p(t, t, prec);
            acb_mul_onei(t, t);
            acb_mul_2exp_si(r, t, -1);
        }
        else
        {
            mag_t err, err2;
            mag_init(err);
            mag_init(err2);

            /* atan'(z) = 1/(1+z^2) = 1/((z+i)(z-i)) */
            acb_onei(t);
            acb_add(t, z, t, prec);
            acb_get_mag_lower(err, t);

            acb_onei(t);
            acb_sub(t, z, t, prec);
            acb_get_mag_lower(err2, t);

            mag_mul_lower(err, err, err2);

            mag_hypot(err2, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
            mag_div(err, err2, err);

            arf_set(arb_midref(acb_realref(t)), arb_midref(acb_realref(z)));
            arf_set(arb_midref(acb_imagref(t)), arb_midref(acb_imagref(z)));
            mag_zero(arb_radref(acb_realref(t)));
            mag_zero(arb_radref(acb_imagref(t)));
            acb_atan(r, t, prec);
            acb_add_error_mag(r, err);

            mag_clear(err);
            mag_clear(err2);
        }

        acb_clear(t);
        acb_clear(u);
    }
}

