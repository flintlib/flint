/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sqrt_ui(arb_t z, ulong x, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_sqrt_arf(z, t, prec);
}

void
arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_sqrt_arf(z, t, prec);
    arf_clear(t);
}

void
arb_sqrt_arf(arb_t z, const arf_t x, slong prec)
{
    if (arf_sgn(x) < 0 || arf_is_nan(x))
    {
        arb_indeterminate(z);
    }
    else
    {
        int inexact;

        inexact = arf_sqrt(arb_midref(z), x, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
}

void
arb_sqrt(arb_t z, const arb_t x, slong prec)
{
    mag_t rx, zr;
    int inexact;

    if (mag_is_zero(arb_radref(x)))
    {
        arb_sqrt_arf(z, arb_midref(x), prec);
    }
    else if (arf_is_special(arb_midref(x)) ||
              arf_sgn(arb_midref(x)) < 0 || mag_is_inf(arb_radref(x)))
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_sqrt_arf(z, arb_midref(x), prec);
        else
            arb_indeterminate(z);
    }
    else  /* now both mid and rad are non-special values, mid > 0 */
    {
        slong acc;

        acc = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)), MAG_EXPREF(arb_radref(x)));
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 0)
        {
            arb_indeterminate(z);
        }
        else if (acc <= 20)
        {
            mag_t t, u;

            mag_init(t);
            mag_init(u);

            arb_get_mag_lower(t, x);

            if (mag_is_zero(t) && arb_contains_negative(x))
            {
                arb_indeterminate(z);
            }
            else
            {
                arb_get_mag(u, x);
                mag_sqrt_lower(t, t);
                mag_sqrt(u, u);
                arb_set_interval_mag(z, t, u, prec);
            }

            mag_clear(t);
            mag_clear(u);
        }
        else if (ARB_IS_LAGOM(x)) /* small exponents, acc *and* prec >= 20 */
        {
            mag_t t;
            mag_init(t); /* no need to free */

            inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

            /* sqrt(x) - sqrt(x-r) <= 0.5 * r * rsqrt(x-r)  */
            /* we have rsqrt(x-r) ~= 1/sqrt(x) */
            arf_get_mag_lower(t, arb_midref(z));

            /* note: we need to write rad(z) first to use fast_mul later */
            mag_div(arb_radref(z), arb_radref(x), t);

            /* We are guaranteed to have acc and prec >= 20. */
            /* 0.5 + eps corrects for errors */
            MAG_MAN(t) = MAG_ONE_HALF + (MAG_ONE_HALF >> 16);
            MAG_EXP(t) = 0;
            mag_fast_mul(arb_radref(z), arb_radref(z), t);

            if (inexact)
                arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        }
        else
        {
            mag_init(zr);
            mag_init(rx);

            /* rx = upper bound for r / x */
            arf_get_mag_lower(rx, arb_midref(x));
            mag_div(rx, arb_radref(x), rx);

            inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

            /* zr = upper bound for sqrt(x) */
            arf_get_mag(zr, arb_midref(z));
            if (inexact)
                arf_mag_add_ulp(zr, zr, arb_midref(z), prec);

            /* propagated error:   sqrt(x) - sqrt(x-r)
                                 = sqrt(x) * [1 - sqrt(1 - r/x)]
                                <= sqrt(x) * 0.5 * (rx + rx^2)  */
            mag_addmul(rx, rx, rx);
            mag_mul(zr, zr, rx);
            mag_mul_2exp_si(zr, zr, -1);

            /* add the rounding error */
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), zr);

            mag_clear(zr);
            mag_clear(rx);
        }
    }
}

