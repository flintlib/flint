/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "double_extras.h"
#include "arf.h"

/* most double: (2^53-1) * 2^971 */
/* least normal: 2^-1022 */
/* least subnormal: 2^-1074 */

static double
huge_double(arf_rnd_t rnd, int negative)
{
    double v;

    if (rnd == ARF_RND_NEAR || arf_rounds_up(rnd, negative))
        v = D_INF;
    else
        v = ldexp(9007199254740991.0, 971);

    return negative ? -v : v;
}

static double
tiny_double(arf_rnd_t rnd, int negative)
{
    double v;

    if (rnd == ARF_RND_NEAR || !arf_rounds_up(rnd, negative))
        v = 0.0;
    else
        v = ldexp(1.0, -1074);

    return negative ? -v : v;
}

double
arf_get_d(const arf_t x, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return 0.0;
        else if (arf_is_pos_inf(x))
            return D_INF;
        else if (arf_is_neg_inf(x))
            return -D_INF;
        else
            return D_NAN;
    }
    else
    {
        arf_t t;
        mp_srcptr tp;
        mp_size_t tn;
        double v;

        /* also catches bignum exponents */
        if (ARF_EXP(x) > 1030 || ARF_EXP(x) < -1080)
        {
            if (fmpz_sgn(ARF_EXPREF(x)) > 0)
                return huge_double(rnd, ARF_SGNBIT(x));
            else
                return tiny_double(rnd, ARF_SGNBIT(x));
        }

        /* allow mpfr to take care of corner cases for now */
        if (ARF_EXP(x) > 1020 || ARF_EXP(x) <= -1020 || rnd == ARF_RND_NEAR)
        {
            mpfr_t xx;
            ARF_GET_MPN_READONLY(tp, tn, x);

            xx->_mpfr_d = (mp_ptr) tp;
            xx->_mpfr_prec = tn * FLINT_BITS;
            xx->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
            xx->_mpfr_exp = ARF_EXP(x);

            return mpfr_get_d(xx, arf_rnd_to_mpfr(rnd));
        }

        arf_init(t);
        arf_set_round(t, x, 53, rnd);
        ARF_GET_MPN_READONLY(tp, tn, t);

        if (tn == 1)
            v = (double)(tp[0]);
        else
            v = (double)(tp[1]) + (double)(tp[0]) * ldexp(1,-32);

        v = ldexp(v, ARF_EXP(t) - FLINT_BITS);

        if (ARF_SGNBIT(t))
            v = -v;

        arf_clear(t);

        return v;
    }
}

