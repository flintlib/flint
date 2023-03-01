/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
arf_shallow_set_siui(arf_t res, ulong vhi, ulong vlo)
{
    int negative;
    unsigned int bc;

    negative = ((slong) vhi) < 0;

    if (negative)
    {
        vhi = -vhi - (vlo != 0);
        vlo = -vlo;
    }

    if (vhi == 0)
    {
        if (vlo == 0)
        {
            ARF_XSIZE(res) = 0;
            ARF_EXP(res) = ARF_EXP_ZERO;
        }
        else
        {
            count_leading_zeros(bc, vlo);
            ARF_EXP(res) = FLINT_BITS - bc;
            ARF_NOPTR_D(res)[0] = vlo << bc;
            ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
        }
    }
    else if (vlo == 0)
    {
        count_leading_zeros(bc, vhi);
        ARF_EXP(res) = 2 * FLINT_BITS - bc;
        ARF_NOPTR_D(res)[0] = vhi << bc;
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
    }
    else
    {
        count_leading_zeros(bc, vhi);
        ARF_EXP(res) = 2 * FLINT_BITS - bc;
        ARF_NOPTR_D(res)[0] = vlo << bc;
        if (bc == 0)
            ARF_NOPTR_D(res)[1] = vhi;
        else
            ARF_NOPTR_D(res)[1] = (vhi << bc) | (vlo >> (FLINT_BITS - bc));
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, negative);
    }
}

void
acb_dot_siui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
{
    arb_ptr t;
    slong i;
    ulong vhi, vlo;
    TMP_INIT;

    /* todo: fast fma and fmma (len=2) code */
    if (len <= 1)
    {
        if (initial == NULL)
        {
            if (len <= 0)
                acb_zero(res);
            else
            {
                arf_t t;
                arf_shallow_set_siui(t, y[1], y[0]);
                arb_mul_arf(acb_realref(res), acb_realref(x), t, prec);
                arb_mul_arf(acb_imagref(res), acb_imagref(x), t, prec);
                if (subtract)
                    acb_neg(res, res);
            }
            return;
        }
        else if (len <= 0)
        {
            acb_set_round(res, initial, prec);
            return;
        }
    }

    TMP_START;
    t = TMP_ALLOC(sizeof(arb_struct) * len);

    for (i = 0; i < len; i++)
    {
        vlo = y[2 * i * ystep];
        vhi = y[2 * i * ystep + 1];

        arf_shallow_set_siui(arb_midref(t + i), vhi, vlo);

        MAG_EXP(arb_radref(t + i)) = 0;
        MAG_MAN(arb_radref(t + i)) = 0;
    }

    arb_dot(((arb_ptr) res) + 0, (initial == NULL) ? NULL : ((arb_srcptr) initial) + 0, subtract, ((arb_srcptr) x) + 0, 2 * xstep, t, 1, len, prec);
    arb_dot(((arb_ptr) res) + 1, (initial == NULL) ? NULL : ((arb_srcptr) initial) + 1, subtract, ((arb_srcptr) x) + 1, 2 * xstep, t, 1, len, prec);

    TMP_END;
}

