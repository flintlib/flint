/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_sub_special(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_zero(x))
    {
        if (arf_is_zero(y))
        {
            arf_zero(z);
            return 0;
        }
        else
            return arf_neg_round(z, y, prec, rnd);
    }
    else if (arf_is_zero(y))
    {
        return arf_set_round(z, x, prec, rnd);
    }
    else if (arf_is_nan(x) || arf_is_nan(y)
        || (arf_is_pos_inf(x) && arf_is_pos_inf(y))
        || (arf_is_neg_inf(x) && arf_is_neg_inf(y)))
    {
        arf_nan(z);
        return 0;
    }
    else if (arf_is_special(x))
    {
        arf_set(z, x);
        return 0;
    }
    else
    {
        arf_neg(z, y);
        return 0;
    }
}

int
arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    mp_srcptr xptr, yptr;
    slong shift;

    if (arf_is_special(x) || arf_is_special(y))
    {
        return arf_sub_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(ARF_EXPREF(x), ARF_EXPREF(y));

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    if (shift >= 0)
        return _arf_add_mpn(z, xptr, xn, ARF_SGNBIT(x), ARF_EXPREF(x),
                           yptr, yn, ARF_SGNBIT(y) ^ 1, shift, prec, rnd);
    else
        return _arf_add_mpn(z, yptr, yn, ARF_SGNBIT(y) ^ 1, ARF_EXPREF(y),
                           xptr, xn, ARF_SGNBIT(x), -shift, prec, rnd);
}

int
arf_sub_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    mp_srcptr xptr, yptr;
    mp_limb_t ytmp;
    int xsgnbit, ysgnbit;
    fmpz yexp;
    slong shift;

    if (y == 0)
    {
        return arf_set_round(z, x, prec, rnd);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arf_set_si(z, y);
            return arf_neg_round(z, z, prec, rnd);
        }
        else
        {
            arf_set(z, x);
            return 0;
        }
    }

    ysgnbit = (y < 0);
    if (ysgnbit)
        ytmp = -y;
    else
        ytmp = y;
    yptr = &ytmp;
    yn = 1;
    yexp = FLINT_BITS;
    ysgnbit ^= 1;

    shift = _fmpz_sub_small(ARF_EXPREF(x), &yexp);

    xsgnbit = ARF_SGNBIT(x);
    ARF_GET_MPN_READONLY(xptr, xn, x);

    if (shift >= 0)
        return _arf_add_mpn(z, xptr, xn, xsgnbit, ARF_EXPREF(x),
                               yptr, yn, ysgnbit, shift, prec, rnd);
    else
        return _arf_add_mpn(z, yptr, yn, ysgnbit, &yexp,
                               xptr, xn, xsgnbit, -shift, prec, rnd);
}

int
arf_sub_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    mp_srcptr xptr, yptr;
    mp_limb_t ytmp;
    int xsgnbit, ysgnbit;
    fmpz yexp;
    slong shift;

    if (y == 0)
    {
        return arf_set_round(z, x, prec, rnd);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arf_set_ui(z, y);
            return arf_neg_round(z, z, prec, rnd);
        }
        else
        {
            arf_set(z, x);
            return 0;
        }
    }

    ysgnbit = 1;
    ytmp = y;
    yptr = &ytmp;
    yn = 1;

    yexp = FLINT_BITS;
    shift = _fmpz_sub_small(ARF_EXPREF(x), &yexp);

    xsgnbit = ARF_SGNBIT(x);
    ARF_GET_MPN_READONLY(xptr, xn, x);

    if (shift >= 0)
        return _arf_add_mpn(z, xptr, xn, xsgnbit, ARF_EXPREF(x),
                               yptr, yn, ysgnbit, shift, prec, rnd);
    else
        return _arf_add_mpn(z, yptr, yn, ysgnbit, &yexp,
                               xptr, xn, xsgnbit, -shift, prec, rnd);
}

int
arf_sub_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    mp_srcptr xptr, yptr;
    mp_limb_t ytmp;
    int xsgnbit, ysgnbit;
    fmpz yexp;
    slong shift;

    if (fmpz_is_zero(y))
    {
        return arf_set_round(z, x, prec, rnd);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arf_set_fmpz(z, y);
            return arf_neg_round(z, z, prec, rnd);
        }
        else
        {
            arf_set(z, x);
            return 0;
        }
    }

    FMPZ_GET_MPN_READONLY(ysgnbit, yn, yptr, ytmp, *y)
    yexp = yn * FLINT_BITS;
    shift = _fmpz_sub_small(ARF_EXPREF(x), &yexp);
    ysgnbit ^= 1;

    xsgnbit = ARF_SGNBIT(x);
    ARF_GET_MPN_READONLY(xptr, xn, x);

    if (shift >= 0)
        return _arf_add_mpn(z, xptr, xn, xsgnbit, ARF_EXPREF(x),
                               yptr, yn, ysgnbit, shift, prec, rnd);
    else
        return _arf_add_mpn(z, yptr, yn, ysgnbit, &yexp,
                               xptr, xn, xsgnbit, -shift, prec, rnd);
}
