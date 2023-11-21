/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define BINEXP_LIMIT 64

typedef enum {POSITIVE = 0, NEGATIVE_EVEN, NEGATIVE_ODD} sign_type;

void
_arb_pow_exp(arb_t z, const arb_t x, sign_type negx, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);

    if (negx == POSITIVE)
        arb_log(t, x, prec);
    else
    {
        arb_neg(t, x);
        arb_log(t, t, prec);
    }

    arb_mul(t, t, y, prec);
    arb_exp(z, t, prec);

    if (negx == NEGATIVE_ODD)
        arb_neg(z, z);

    arb_clear(t);
}

void
arb_pow(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    sign_type s;

    if (arb_is_zero(y))
    {
        arb_one(z);
        return;
    }

    if (arb_is_zero(x))
    {
        if (arb_is_positive(y))
            arb_zero(z);
        else
            arb_indeterminate(z);
        return;
    }


    s = POSITIVE;

    if (arb_is_exact(y) && !arf_is_special(arb_midref(x)))
    {
        const arf_struct * ymid = arb_midref(y);

        /* small half-integer or integer */
        if (arf_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0 &&
            arf_is_int_2exp_si(ymid, -1))
        {
            fmpz_t e;
            fmpz_init(e);

            if (arf_is_int(ymid))
            {
                arf_get_fmpz_fixed_si(e, ymid, 0);
                arb_pow_fmpz_binexp(z, x, e, prec);
            }
            else
            {
                arf_get_fmpz_fixed_si(e, ymid, -1);
                if (fmpz_sgn(e) >= 0)
                {
                    arb_sqrt(z, x, prec + fmpz_bits(e));
                    arb_pow_fmpz_binexp(z, z, e, prec);
                }
                else
                {
                    fmpz_neg(e, e);
                    arb_rsqrt(z, x, prec + fmpz_bits(e));
                    arb_pow_fmpz_binexp(z, z, e, prec);
                }
            }

            fmpz_clear(e);
            return;
        }
        else if (arf_is_int(ymid) && arf_sgn(arb_midref(x)) < 0)
            /* use (-x)^n = (-1)^n * x^n to avoid NaNs
               at least at high enough precision */
            s = arf_is_int_2exp_si(ymid, 1) ? NEGATIVE_EVEN : NEGATIVE_ODD;
        /* Fallthrough */
    }

    if (arf_cmp_si(arb_midref(x), 0) > 0
        && arf_cmpabs_mag(arb_midref(x), arb_radref(x)) == 0
        && arb_is_nonnegative(y)) {
        /* x is an interval of the form <a +- a>, y is an interval of the form
           <b +- c> with b >= c.

           (x, y) -> x^y is nondecreasing in x for y >= 0, and it is 0 for x=0
           and y>0.  The lower bound of the interval for x is 0, and we have
           points with y>0 (because the case arb_is_zero(y) is excluded
           above). So the lowerbound of the result is 0, and the upperbound can
           be found somewhere at (2*a)^y.

           We compute this upperbound by computing z^y with z = <3*a/2 +- a/2>,
           which has the same upperbound (2*a) as x, then keeping the upper
           bound of this result and setting the lower bound to 0. This
           necessarily drops the precision to MAG_BITS, so we might as well do
           that from the start. */
        prec = (prec > MAG_BITS) ? MAG_BITS : prec;
        arf_mul_ui(arb_midref(z), arb_midref(x), 3, prec, ARF_RND_UP);
        arf_mul_2exp_si(arb_midref(z), arb_midref(z), -1);
        mag_mul_2exp_si(arb_radref(z), arb_radref(x), -1);

        arb_pow(z, z, y, prec);

        /* Now we keep the upper bound of z (we may need to round up) and set the lower bound to
           zero. That is, if currently, z = <zc +/- zr>, we set it to <(zc + zr)/2 +/- (zc +
           zr)/2>. */
        arb_get_ubound_arf(arb_midref(z), z, prec);
        arf_mul_2exp_si(arb_midref(z), arb_midref(z), -1);
        arf_get_mag(arb_radref(z), arb_midref(z));
        /* Note the above is inexact (rounding up), so need to update arb_midref(z) to match again */
        arf_set_mag(arb_midref(z), arb_radref(z));

        return;
    }

    _arb_pow_exp(z, x, s, y, prec);
}

