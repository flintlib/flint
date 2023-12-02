/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpfr.h"
#include "arf.h"

typedef int ((*mpfr_func_1x1)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t));
typedef int ((*mpfr_func_1x2)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t));
typedef int ((*mpfr_func_2x1)(mpfr_ptr, mpfr_ptr, mpfr_srcptr, mpfr_rnd_t));

int _arf_call_mpfr_func(arf_ptr r1, arf_ptr r2, int (*func)(void), arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    mpfr_t xx, yy, rr1, rr2;
    mpfr_rnd_t rrnd;
    int inexact = 0;

    rrnd = arf_rnd_to_mpfr(rnd);

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    mpfr_init2(xx, 2 + arf_bits(x));
    if (arf_get_mpfr(xx, x, MPFR_RNDD))
    {
        flint_throw(FLINT_ERROR, "exception: unable to convert exactly to mpfr\n");
    }

    if (y != NULL)
    {
        mpfr_init2(yy, 2 + arf_bits(y));
        if (arf_get_mpfr(yy, y, MPFR_RNDD))
        {
            flint_throw(FLINT_ERROR, "exception: unable to convert exactly to mpfr\n");
        }
    }

    mpfr_init2(rr1, FLINT_MAX(2, prec));

    if (r2 != NULL)
        mpfr_init2(rr2, FLINT_MAX(2, prec));

    if (r2 == NULL && y == NULL)
        inexact = (((mpfr_func_1x1) func)(rr1, xx, rrnd) != 0);
    else if (r2 != NULL && y == NULL)
        inexact = (((mpfr_func_2x1) func)(rr1, rr2, xx, rrnd) != 0);
    else if (r2 == NULL && y != NULL)
        inexact = (((mpfr_func_1x2) func)(rr1, xx, yy, rrnd) != 0);
    else
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (mpfr_overflow_p() || mpfr_underflow_p())
    {
        flint_throw(FLINT_ERROR, "exception: mpfr overflow\n");
    }

    if (r1 != NULL)
    {
        arf_set_mpfr(r1, rr1);
        mpfr_clear(rr1);
    }

    if (r2 != NULL)
    {
        arf_set_mpfr(r2, rr2);
        mpfr_clear(rr2);
    }

    if (x != NULL)
        mpfr_clear(xx);
    if (y != NULL)
        mpfr_clear(yy);

    return inexact;
}
