/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

void
_mp_real_elem_half_ratio(nn_ptr w, nn_srcptr s, slong n, int minus)
{
    nn_ptr D;
    TMP_INIT;

    if (n <= 2)
    {
        /* the exact floor, in registers */
        ulong s1 = (n == 2) ? s[1] : s[0], s0 = (n == 2) ? s[0] : 0;
        ulong d2, d1, d0, q[2];

        if (minus)
        {
            sub_dddmmmsss(d2, d1, d0, UWORD(2), UWORD(0), UWORD(0), UWORD(0), s1, s0);
        }
        else
        {
            d2 = 2; d1 = s1; d0 = s0;
        }

        _mp_real_divq_4_3z(q, s1, s0, d2, d1, d0);

        if (n == 2)
        {
            w[0] = q[0];
            w[1] = q[1];
        }
        else
            w[0] = q[1];
        return;
    }

    TMP_START;
    D = TMP_ALLOC((n + 1) * sizeof(ulong));
    if (minus)
    {
        /* 2 B^n - s = B^n + (B^n - s) */
        mpn_neg(D, s, n);
        D[n] = 1;
    }
    else
    {
        flint_mpn_copyi(D, s, n);
        D[n] = 2;
    }
    /* floor(s B^n / D) or one more: n limbs, below B^n since s < D */
    flint_mpn_divapprox_fraction(w, s, n, D, n + 1, n);
    TMP_END;
}

/* Normalisation of fixed-point results *************************************/

ulong
_mp_real_fp_normalize(nn_ptr t, slong * te, nn_srcptr y, slong len, slong ue, ulong err, slong wn)
{
    slong sh, i;
    unsigned int c;

    while (y[len - 1] == 0)
        len--;

    c = flint_clz(y[len - 1]);
    sh = FLINT_BITS * len - c - FLINT_BITS * wn;    /* right shift */
    *te = ue + sh;

    if (sh >= 0)
    {
        slong ls = sh / FLINT_BITS;
        unsigned int bs = sh % FLINT_BITS;
        int dropped = 0;

        for (i = 0; i < ls; i++)
            dropped |= (y[i] != 0);

        if (bs == 0)
            _flint_mpn_copy_small(t, y + ls, wn);
        else
        {
            dropped |= ((y[ls] << (FLINT_BITS - bs)) != 0);
            for (i = 0; i < wn; i++)
                t[i] = (y[ls + i] >> bs) | ((ls + i + 1 < len) ? (y[ls + i + 1] << (FLINT_BITS - bs)) : 0);
        }

        if (sh >= FLINT_BITS)
            err = 1;
        else
            err = (err >> sh) + (sh != 0);

        return err + dropped;
    }
    else
    {
        slong ls = (-sh) / FLINT_BITS;
        unsigned int bs = (-sh) % FLINT_BITS;

        _flint_mpn_zero_small(t, ls);
        if (bs == 0)
            _flint_mpn_copy_small(t + ls, y, len);
        else
        {
            t[ls] = y[0] << bs;
            for (i = 1; i < len; i++)
                t[ls + i] = (y[i] << bs) | (y[i - 1] >> (FLINT_BITS - bs));
        }

        /* the error scales up by 2^-sh (UWORD_MAX signals overflow) */
        if (-sh >= FLINT_BITS - 2 || err > (UWORD(1) << (FLINT_BITS - 3 + sh)))
            return UWORD_MAX;
        return err << (-sh);
    }
}
