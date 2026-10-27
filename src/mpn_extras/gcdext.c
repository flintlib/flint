/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

void
_flint_mpn_gcdext_cofactor(mp_ptr tp, mp_size_t * tn,
    mp_srcptr gp, mp_size_t gn, mp_srcptr sp, mp_size_t sn,
    mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_srcptr xp;
    mp_ptr yp, qp;
    mp_size_t xn, qn, abs_sn = FLINT_ABS(sn);
    int xneg;
    TMP_INIT;

    TMP_START;

    if (sn == 0)
    {
        /* x = g */
        xp = gp;
        xn = gn;
        xneg = 0;
    }
    else
    {
        yp = TMP_ALLOC((an + abs_sn + 1) * sizeof(mp_limb_t));
        if (an >= abs_sn)
            flint_mpn_mul(yp, ap, an, sp, abs_sn);
        else
            flint_mpn_mul(yp, sp, abs_sn, ap, an);
        xn = an + abs_sn;
        xn -= (yp[xn - 1] == 0);

        if (sn < 0)
        {
            /* x = g + |s| a > 0 */
            mp_limb_t cy = mpn_add(yp, yp, xn, gp, gn);
            yp[xn] = cy;
            xn += (cy != 0);
            xneg = 0;
        }
        else
        {
            /* x = g - s a <= 0; |x| = s a - g */
            mpn_sub(yp, yp, xn, gp, gn);
            MPN_NORM(yp, xn);
            xneg = 1;
        }

        xp = yp;
    }

    if (xn == 0)
    {
        *tn = 0;
    }
    else
    {
        FLINT_ASSERT(xn >= bn);
        qn = xn - bn + 1;
        qp = TMP_ALLOC(qn * sizeof(mp_limb_t));
        flint_mpn_divexact(qp, xp, xn, bp, bn);
        MPN_NORM(qp, qn);
        FLINT_ASSERT(qn <= an);
        flint_mpn_copyi(tp, qp, qn);
        *tn = xneg ? -qn : qn;
    }

    TMP_END;
}

mp_size_t
flint_mpn_gcdext(mp_ptr gp, mp_ptr sp, mp_size_t * sn,
    mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_ptr tp, ta, tb, tg, ts;
    mp_size_t gn, tsn, m, n;
    TMP_INIT;

    FLINT_ASSERT(an > 0 && ap[an - 1] != 0);
    FLINT_ASSERT(bn > 0 && bp[bn - 1] != 0);

    /* the inplace version requires that its first operand is at least as
       long as the second; it then returns the cofactor of the first */
    m = FLINT_MAX(an, bn);
    n = FLINT_MIN(an, bn);

    TMP_START;

    tp = TMP_ALLOC((m + 1 + n + n + n + 1) * sizeof(mp_limb_t));
    ta = tp;
    tb = ta + m + 1;
    tg = tb + n;
    ts = tg + n;

    if (an >= bn)
    {
        flint_mpn_copyi(ta, ap, an);
        flint_mpn_copyi(tb, bp, bn);
        gn = _flint_mpn_gcdext_inplace(tg, ts, &tsn, ta, an, tb, bn);
        flint_mpn_copyi(sp, ts, FLINT_ABS(tsn));
        *sn = tsn;
    }
    else
    {
        /* the cofactor of b, then that of a = (g - s_b b) / a */
        flint_mpn_copyi(ta, bp, bn);
        flint_mpn_copyi(tb, ap, an);
        gn = _flint_mpn_gcdext_inplace(tg, ts, &tsn, ta, bn, tb, an);
        _flint_mpn_gcdext_cofactor(sp, sn, tg, gn, ts, tsn, bp, bn, ap, an);
    }

    /* written last, as gp may be ap or bp */
    flint_mpn_copyi(gp, tg, gn);

    TMP_END;
    return gn;
}
