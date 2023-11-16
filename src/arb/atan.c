/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_atan(arb_t res, const arb_t x, slong prec)
{
    if (mag_is_zero(arb_radref(x)))
    {
        arb_atan_arf(res, arb_midref(x), prec);
    }
    else if (arf_is_nan(arb_midref(x)))
    {
        arb_indeterminate(res);
    }
    else if (mag_is_inf(arb_radref(x)) || arf_is_zero(arb_midref(x)))
    {
        mag_atan(arb_radref(res), arb_radref(x));
        arf_zero(arb_midref(res));
    }
    else if (arf_is_special(arb_midref(x)))  /* at +/- inf */
    {
        arb_atan_arf(res, arb_midref(x), prec);
    }
    else  /* both mid(x), rad(x) non-special */
    {
        slong acc;
        mag_t t, u;

        /* Estimate accuracy (first guess, only good near 0). */
        acc = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)),
                              MAG_EXPREF(arb_radref(x)));

        if (acc < -10)   /* Only compute a rough bound. */
        {
            arb_get_mag(arb_radref(res), x);
            mag_atan(arb_radref(res), arb_radref(res));
            arf_zero(arb_midref(res));
            return;
        }

        mag_init(t);
        mag_init(u);

        arb_get_mag_lower(t, x);

        if (mag_is_zero(t))   /* Interval includes zero. */
        {
            /* atan(rad(x) - |mid(x)|) */
            arf_get_mag_lower(t, arb_midref(x));
            mag_sub(t, arb_radref(x), t);
            mag_atan(t, t);

            /* atan(|mid(x)| + rad(x)) */
            arf_get_mag(u, arb_midref(x));
            mag_add(u, arb_radref(x), u);
            mag_atan(u, u);

            if (arf_sgn(arb_midref(x)) > 0)
                arb_set_interval_neg_pos_mag(res, t, u, FLINT_MIN(prec, MAG_BITS));
            else
                arb_set_interval_neg_pos_mag(res, u, t, FLINT_MIN(prec, MAG_BITS));
        }
        else
        {
            /* Adjust estimate of accuracy for large x. */
            if (fmpz_sgn(MAG_EXPREF(t)) > 0)
            {
                /* acc = 2 * texp - radexp */
                acc = _fmpz_sub_small(MAG_EXPREF(t), MAG_EXPREF(arb_radref(x)));
                if (acc < prec && !COEFF_IS_MPZ(MAG_EXP(t)))
                    acc += MAG_EXP(t);
            }

            /* Clamp acc and adjust precision. */
            acc = FLINT_MAX(acc, 0);
            acc = FLINT_MIN(acc, prec);
            prec = FLINT_MIN(prec, acc + MAG_BITS);
            prec = FLINT_MAX(prec, 2);

            /* Compute using endpoints */
            if (acc < 20)
            {
                arb_get_mag(u, x);
                mag_atan_lower(t, t);
                mag_atan(u, u);

                if (arf_sgn(arb_midref(x)) > 0)
                {
                    arb_set_interval_mag(res, t, u, prec);
                }
                else
                {
                    arb_set_interval_mag(res, t, u, prec);
                    arb_neg(res, res);
                }
            }
            else
            {
                mag_mul_lower(t, t, t);
                mag_one(u);
                mag_add_lower(t, t, u);
                mag_div(t, arb_radref(x), t);

                if (mag_cmp_2exp_si(t, 0) > 0)
                {
                    mag_const_pi(u);
                    mag_min(t, t, u);
                }

                arb_atan_arf(res, arb_midref(x), prec);
                mag_add(arb_radref(res), arb_radref(res), t);
            }
        }

        mag_clear(t);
        mag_clear(u);
    }
}

