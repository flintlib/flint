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
arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, zn, val;
    mp_srcptr xptr;
    mp_ptr tmp, zptr;
    mpfr_t xf, zf;
    int inexact, odd_exp;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_pos_inf(z);
        else if (arf_is_pos_inf(x))
            arf_zero(z);
        else
            arf_nan(z);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_set(z, x);
        fmpz_neg(ARF_EXPREF(z), ARF_EXPREF(z));
        fmpz_cdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        fmpz_add_ui(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);
    odd_exp = fmpz_is_odd(ARF_EXPREF(x)) ? 1 : 0;

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ARF_MUL_TMP_ALLOC(tmp, zn)

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = odd_exp;

    inexact = mpfr_rec_sqrt(zf, xf, arf_rnd_to_mpfr(rnd));
    inexact = (inexact != 0);

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);

    fmpz_fdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(x), 1);
    fmpz_neg(ARF_EXPREF(z), ARF_EXPREF(z));

    fmpz_add_si(ARF_EXPREF(z), ARF_EXPREF(z), zf->_mpfr_exp);

    ARF_MUL_TMP_FREE(tmp, zn)

    return inexact;
}

