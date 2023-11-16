/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define ldexp __builtin_ldexp
# define log __builtin_log
#else
# include <math.h>
#endif

static double
_arf_get_mantissa_d(const arf_t x)
{
    mp_srcptr xp;
    mp_size_t xn;
    ARF_GET_MPN_READONLY(xp, xn, x);

    if (xn == 1)
        return (double) xp[0] * ldexp(1.0, -FLINT_BITS);
    else
        return ((double) xp[xn - 1]) * ldexp(1.0, -FLINT_BITS) +
               ((double) xp[xn - 2]) * ldexp(1.0, -2 * FLINT_BITS);
}

void
arb_log_base_ui(arb_t res, const arb_t x, ulong b, slong prec)
{
    arb_t t;
    slong xexp, xbits;

    if (b <= 1)
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_is_exact(x) && arf_sgn(arb_midref(x)) > 0)
    {
        xbits = arb_bits(x);
        xexp = ARF_EXP(arb_midref(x));

        /* x = 1 */
        if (xbits == 1 && ARF_EXP(arb_midref(x)) == 1)
        {
            arb_zero(res);
            return;
        }

        /* powers of two */
        if ((b & (b - 1)) == 0 && xbits == 1)
        {
            fmpz_t e;
            fmpz_init(e);
            fmpz_sub_ui(e, ARF_EXPREF(arb_midref(x)), 1);

            if (b == 2)
            {
                arb_set_round_fmpz(res, e, prec);
            }
            else
            {
                arb_set_fmpz(res, e);
                arb_div_ui(res, res, FLINT_BIT_COUNT(b) - 1, prec);
            }
            fmpz_clear(e);
            return;
        }

        /* check for other exact powers */
        if ((b & (b - 1)) != 0 && xbits != 1)
        {
            /* exactness is only possible if x is moderate, and an integer */
            if (!COEFF_IS_MPZ(xexp) && xbits <= xexp)
            {
                ulong b_reduced;
                int b_exp;

                /* Write b as (b_reduced)^(1/2^(b_exp)). */
                /* More generally, we could look for nth roots, but the result
                   will not be exact when dividing by a non-power-of-two,
                   so that is just a possible optimization at high precision. */
                b_reduced = b;
                b_exp = 0;

                if (b >= 25 || b == 9)
                {
                    ulong s, r;
                    s = n_sqrtrem(&r, b_reduced);
                    while (r == 0)
                    {
                        b_reduced = s;
                        b_exp++;
                        s = n_sqrtrem(&r, b_reduced);
                    }
                }

                /* x fits a ulong */
                if (xexp <= FLINT_BITS)
                {
                    ulong n, a, hi, v;

                    ARF_GET_TOP_LIMB(v, arb_midref(x));
                    v >>= (FLINT_BITS - xexp);

                    for (a = b_reduced, n = 1, hi = 0; a <= v && hi == 0; n++)
                    {
                        if (a == v)
                        {
                            arf_set_ui_2exp_si(arb_midref(res), n, -b_exp);
                            mag_zero(arb_radref(res));
                            arb_set_round(res, res, prec);
                            return;
                        }

                        umul_ppmm(hi, a, a, b_reduced);
                    }
                }
                /* general case (if b = 10, first count bits to rule out many
                   non-powers-of-10 -- could perhaps generalize this to other b) */
                else if (b_reduced != 10 || (xbits < xexp && xexp < xbits * 2))
                {
                    ulong n;
                    double xlog;

                    /* libm.log should be accurate enough since we will certainly
                       have bits(n) << 53. Worst case, we will not get the exact
                       logarithm and fall back to compute an approximate logarithm. */
                    xlog = _arf_get_mantissa_d(arb_midref(x));
                    xlog = log(xlog) + xexp * 0.69314718055994530942;
                    xlog = xlog / log(b_reduced);
                    n = (ulong) (xlog + 0.5);

                    if (n >= 2 && fabs(xlog - n) < 0.01)
                    {
                        arb_init(t);
                        arb_ui_pow_ui(t, b_reduced, n, xbits + 10);

                        if (arb_equal(t, x))
                        {
                            arf_set_ui_2exp_si(arb_midref(res), n, -b_exp);
                            mag_zero(arb_radref(res));
                            arb_set_round(res, res, prec);
                            arb_clear(t);
                            return;
                        }

                        arb_clear(t);
                    }
                }
            }
        }
    }

    /* generic case */
    arb_init(t);
    arb_log(res, x, prec + 3);
    arb_log_ui(t, b, prec + 3);
    arb_div(res, res, t, prec);
    arb_clear(t);
}

