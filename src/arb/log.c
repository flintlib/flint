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
arb_log(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_exact(x))
    {
        arb_log_arf(res, arb_midref(x), prec);
    }
    else
    {
        /*
        Let the input be [a-b, a+b]. We require a > b >= 0 (otherwise the
        interval contains zero or a negative number and the logarithm is not
        defined). The error is largest at a-b, and we have

        log(a) - log(a-b) = log(1 + b/(a-b)).
        */
        mag_t t;
        mag_init(t);

        arb_get_mag_lower_nonnegative(t, x);

        if (mag_is_zero(t))
        {
            arf_nan(arb_midref(res));
            mag_inf(arb_radref(res));
        }
        else if (mag_is_inf(t))
        {
            arf_pos_inf(arb_midref(res));
            mag_zero(arb_radref(res));
        }
        else
        {
            slong acc;

            acc = arb_rel_accuracy_bits(x);
            acc = FLINT_MIN(acc, prec);
            acc += fmpz_bits(MAG_EXPREF(t));

            if (acc < 20)
            {
                mag_t u;
                mag_init(u);
                arb_get_mag(u, x);

                if (mag_cmp_2exp_si(t, 0) >= 0)
                {
                    mag_log_lower(t, t);
                    mag_log(u, u);
                    arb_set_interval_mag(res, t, u, prec);
                }
                else if (mag_cmp_2exp_si(u, 0) <= 0)
                {
                    mag_neg_log_lower(u, u);
                    mag_neg_log(t, t);
                    arb_set_interval_mag(res, u, t, prec);
                    arb_neg(res, res);
                }
                else
                {
                    mag_neg_log(t, t);
                    mag_log(u, u);
                    arb_set_interval_neg_pos_mag(res, t, u, prec);
                }

                mag_clear(u);
            }
            else
            {
                acc = FLINT_MAX(acc, 0);
                acc = FLINT_MIN(acc, prec);
                prec = FLINT_MIN(prec, acc + MAG_BITS);

                mag_div(t, arb_radref(x), t);
                mag_log1p(t, t);
                arb_log_arf(res, arb_midref(x), prec);
                mag_add(arb_radref(res), arb_radref(res), t);
            }
        }

        mag_clear(t);
    }
}
