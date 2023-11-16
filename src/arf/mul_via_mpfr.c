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
arf_mul_via_mpfr(arf_t z, const arf_t x, const arf_t y,
        slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, val;
    mp_srcptr xptr, yptr;
    mp_ptr tmp, zptr;
    mpfr_t xf, yf, zf;
    int ret;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y))
    {
        arf_mul_special(z, x, y);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    prec = FLINT_MIN((xn + yn) * FLINT_BITS, prec);
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ARF_MUL_TMP_ALLOC(tmp, zn)

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
    xf->_mpfr_exp = 0;

    if (x == y)
    {
        ret = mpfr_sqr(zf, xf, arf_rnd_to_mpfr(rnd));
    }
    else
    {
        yf->_mpfr_d = (mp_ptr) yptr;
        yf->_mpfr_prec = yn * FLINT_BITS;
        yf->_mpfr_sign = ARF_SGNBIT(y) ? -1 : 1;
        yf->_mpfr_exp = 0;

        ret = mpfr_mul(zf, xf, yf, arf_rnd_to_mpfr(rnd));
    }

    ret = (ret != 0);

    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), zf->_mpfr_exp);

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    ARF_XSIZE(z) |= (zf->_mpfr_sign < 0);

    ARF_MUL_TMP_FREE(tmp, zn)

    return ret;
}

