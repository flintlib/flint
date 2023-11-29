/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arf.h"

int
arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd)
{
    int r;

    if (arf_is_special(y))
    {
        if (arf_is_zero(y))
            mpfr_set_zero(x, 0);
        else if (arf_is_pos_inf(y))
            mpfr_set_inf(x, 1);
        else if (arf_is_neg_inf(y))
            mpfr_set_inf(x, -1);
        else
            mpfr_set_nan(x);
        r = 0;
    }
    else if (COEFF_IS_MPZ(*ARF_EXPREF(y)))
    {
        /* Incidentally, COEFF_MIN and COEFF_MAX are exactly the same
           as the minimum and maximum allowed MPFR exponents. We
           assert this to make sure the following code is valid.
           Unfortunately, MPFR provides no convenient function to
           assign a big exponent with automatic underflow/overflow. */
        if (COEFF_MIN > mpfr_get_emin_min() ||
            COEFF_MAX < mpfr_get_emax_max())
        {
            flint_throw(FLINT_ERROR, "unsupported MPFR exponent range: %wd, %wd, %wd, %wd\n",
                COEFF_MIN, mpfr_get_emin_min(), COEFF_MAX, mpfr_get_emax_max());
        }

        if (fmpz_sgn(ARF_EXPREF(y)) > 0)
        {
            if (arf_sgn(y) > 0)
            {
                mpfr_set_inf(x, 1);
                mpfr_nextbelow(x);
            }
            else
            {
                mpfr_set_inf(x, -1);
                mpfr_nextabove(x);
            }

            r = mpfr_mul_2si(x, x, 1, rnd);
        }
        else
        {
            if (arf_sgn(y) > 0)
            {
                mpfr_set_zero(x, 1);
                mpfr_nextabove(x);
            }
            else
            {
                mpfr_set_zero(x, -1);
                mpfr_nextbelow(x);
            }

            r = mpfr_mul_2si(x, x, -1, rnd);
        }
    }
    else
    {
        __mpfr_struct t;
        mp_size_t n;
        mp_srcptr d;

        ARF_GET_MPN_READONLY(d, n, y);

        t._mpfr_d = (mp_ptr) d;
        t._mpfr_exp = ARF_EXP(y);
        t._mpfr_prec = n * FLINT_BITS;
        t._mpfr_sign = ARF_SGNBIT(y) ? -1 : 1;

        r = mpfr_set(x, &t, rnd);
    }

    return r;
}

