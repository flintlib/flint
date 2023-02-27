/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"
int
fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd)
{
    int r;

    if (fmpr_is_special(y))
    {
        if (fmpr_is_zero(y)) mpfr_set_zero(x, 0);
        else if (fmpr_is_pos_inf(y)) mpfr_set_inf(x, 1);
        else if (fmpr_is_neg_inf(y)) mpfr_set_inf(x, -1);
        else mpfr_set_nan(x);
        r = 0;
    }
    else if (COEFF_IS_MPZ(*fmpr_expref(y)))
    {
        flint_printf("exception: exponent too large to convert to mpfr");
        flint_abort();
        r = 0; /* dummy return because flint_abort() is not declared noreturn */
    }
    else
    {
        if (!COEFF_IS_MPZ(*fmpr_manref(y)))
#if defined(__MINGW64__) || defined(_MSC_VER)
            r = mpfr_set_sj_2exp(x, *fmpr_manref(y), *fmpr_expref(y), rnd);
#else
            r = mpfr_set_si_2exp(x, *fmpr_manref(y), *fmpr_expref(y), rnd);
#endif
        else
            r = mpfr_set_z_2exp(x, COEFF_TO_PTR(*fmpr_manref(y)), *fmpr_expref(y), rnd);

        if (!mpfr_regular_p(x))
        {
            flint_printf("exception: exponent too large to convert to mpfr");
            flint_abort();
        }
    }

    return r;
}

