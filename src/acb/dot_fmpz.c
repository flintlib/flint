/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_dot_fmpz(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec)
{
    arb_ptr t;
    slong i, ssize, size, tmp_size;
    mp_ptr ztmp;
    fmpz v;
    ulong av, al;
    unsigned int bc;
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
                acb_mul_fmpz(res, x, y, prec);
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

    tmp_size = 0;
    for (i = 0; i < len; i++)
    {
        v = y[i * ystep];

        MAG_EXP(arb_radref(t + i)) = 0;
        MAG_MAN(arb_radref(t + i)) = 0;

        if (v == 0)
        {
            ARF_XSIZE(arb_midref(t + i)) = 0;
            ARF_EXP(arb_midref(t + i)) = ARF_EXP_ZERO;
        }
        else if (!COEFF_IS_MPZ(v))
        {
            av = FLINT_ABS(v);
            bc = flint_clz(av);

            ARF_EXP(arb_midref(t + i)) = FLINT_BITS - bc;
            ARF_NOPTR_D(arb_midref(t + i))[0] = av << bc;
            ARF_XSIZE(arb_midref(t + i)) = ARF_MAKE_XSIZE(1, v < 0);
        }
        else
        {
            __mpz_struct * z = COEFF_TO_PTR(v);

            ssize = z->_mp_size;
            size = FLINT_ABS(ssize);

            av = z->_mp_d[size - 1];
            bc = flint_clz(av);

            if (size == 1)
            {
                ARF_EXP(arb_midref(t + i)) = FLINT_BITS - bc;
                ARF_NOPTR_D(arb_midref(t + i))[0] = av << bc;
                ARF_XSIZE(arb_midref(t + i)) = ARF_MAKE_XSIZE(1, ssize < 0);
            }
            else if (size == 2)
            {
                al = z->_mp_d[0];

                ARF_EXP(arb_midref(t + i)) = 2 * FLINT_BITS - bc;

                if (bc != 0)
                {
                    av = (av << bc) | (al >> (FLINT_BITS - bc));
                    al = al << bc;
                }

                ARF_NOPTR_D(arb_midref(t + i))[0] = al;
                ARF_NOPTR_D(arb_midref(t + i))[1] = av;
                ARF_XSIZE(arb_midref(t + i)) = ARF_MAKE_XSIZE(2, ssize < 0);
            }
            else
            {
                if (bc != 0)
                {
                    tmp_size += size;
                    /* use to flag tmp where we need tmp storage */
                    MAG_MAN(arb_radref(t + i)) = bc;
                }

                ARF_EXP(arb_midref(t + i)) = size * FLINT_BITS - bc;
                ARF_PTR_D(arb_midref(t + i)) = z->_mp_d;
                ARF_XSIZE(arb_midref(t + i)) = ARF_MAKE_XSIZE(size, ssize < 0);
            }
        }
    }

    if (tmp_size != 0)
    {
        ztmp = TMP_ALLOC(sizeof(mp_limb_t) * tmp_size);

        for (i = 0; i < len; i++)
        {
            bc = MAG_MAN(arb_radref(t + i));

            if (bc != 0)
            {
                size = ARF_SIZE(arb_midref(t + i));

                mpn_lshift(ztmp, ARF_PTR_D(arb_midref(t + i)), size, bc);
                ARF_PTR_D(arb_midref(t + i)) = ztmp;
                ztmp += size;
            }

            MAG_MAN(arb_radref(t + i)) = 0;
        }
    }

    arb_dot(((arb_ptr) res) + 0, (initial == NULL) ? NULL : ((arb_srcptr) initial) + 0, subtract, ((arb_srcptr) x) + 0, 2 * xstep, t, 1, len, prec);
    arb_dot(((arb_ptr) res) + 1, (initial == NULL) ? NULL : ((arb_srcptr) initial) + 1, subtract, ((arb_srcptr) x) + 1, 2 * xstep, t, 1, len, prec);

    TMP_END;
}

