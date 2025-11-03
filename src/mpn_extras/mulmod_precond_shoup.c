/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* adapted from flint_mpn_divrem_preinvn */
static mp_limb_t flint_mpn_div_preinvn_4_2(mp_ptr qp, mp_srcptr ap,
                          mp_srcptr d, mp_srcptr dinv)
{
    mp_limb_t cy, hi = 0;
    mp_limb_t t[4];
    mp_limb_t r[4];

    if (mpn_cmp(ap + 2, d, 2) >= 0)
    {
        sub_ddmmss(r[3], r[2], ap[3], ap[2], d[1], d[0]);
        hi = 1;
    }
    else
    {
        r[2] = ap[2];
        r[3] = ap[3];
    }

    FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], dinv[1], dinv[0], r[3], r[2]);
    add_sssaaaaaa(cy, qp[1], qp[0], 0, t[3], t[2], 0, r[3], r[2]);
    FLINT_MPN_MUL_3P2X2(t[2], t[1], t[0], d[1], d[0], qp[1], qp[0]);
    sub_dddmmmsss(cy, r[1], r[0], r[2], ap[1], ap[0], t[2], t[1], t[0]);

    while (cy > 0)
    {
        sub_dddmmmsss(cy, r[1], r[0], cy, r[1], r[0], 0, d[1], d[0]);
        add_ssaaaa(qp[1], qp[0], qp[1], qp[0], 0, 1);
    }

    if (mpn_cmp(r, d, 2) >= 0)
    {
        /* We don't need the remainder */
        /* sub_ddmmss(r[1], r[0], r[1], r[0], d[1], d[0]); */
        /* FLINT_ASSERT(mpn_cmp(r, d, n) < 0); */
        add_ssaaaa(qp[1], qp[0], qp[1], qp[0], 0, 1);
    }

    return hi;
}

void
flint_mpn_mulmod_precond_shoup_precompute(mp_ptr apre, mp_srcptr a, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)
{
    if (n == 2)
    {
        mp_limb_t t[4];
        t[0] = 0;
        t[1] = 0;
        if (norm == 0)
        {
            t[2] = a[0];
            t[3] = a[1];
        }
        else
        {
            t[2] = a[0] << norm;
            t[3] = (a[1] << norm) | (a[0] >> (FLINT_BITS - norm));
        }
        flint_mpn_div_preinvn_4_2(apre, t, dnormed, dinv);
    }
    else
    {
        mp_ptr t;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC(sizeof(mp_limb_t) * (2 * n));

        flint_mpn_zero(t, n);
        if (norm == 0)
            mpn_copyi(t + n, a, n);
        else
            mpn_lshift(t + n, a, n, norm);
        flint_mpn_divrem_preinvn(apre, t, t, 2 * n, dnormed, n, dinv);

        TMP_END;
    }
}

/*
   The original version of Shoup's algorithm uses the exact high
   product, for which we would need to inspect the return limb
   and recompute the low part when rounding cannot be guaranteed.
   However, it can be shown that the approximation computed by
   flint_mpn_mulhigh_n suffices when norm != 0.

   When norm == 0, we can either do an exact mulhigh OR do
   two adjustments (warning: the second adjustment may be too rare to be
   reached by test code).
*/
void
flint_mpn_mulmod_precond_shoup(mp_ptr res, mp_srcptr a, mp_srcptr apre, mp_srcptr b, mp_size_t n, mp_srcptr d, ulong norm)
{
    if (n == 2)
    {
        if (norm != 0)
        {
            mp_limb_t t3, t2, t1, t0;

            FLINT_MPN_MUL_2X2(t3, t2, t1, t0, apre[1], apre[0], b[1], b[0]);
            FLINT_MPN_MULLOW_2X2(t1, t0, t3, t2, d[1], d[0]);
            FLINT_MPN_MULLOW_2X2(t3, t2, a[1], a[0], b[1], b[0]);
            sub_ddmmss(t3, t2, t3, t2, t1, t0);
            if (t3 > d[1] || (t3 == d[1] && t2 >= d[0]))
                sub_ddmmss(t3, t2, t3, t2, d[1], d[0]);
            res[0] = t2;
            res[1] = t3;
        }
        else
        {
            mp_limb_t cy, tcy, ucy, t3, t2, t1, t0;

            FLINT_MPN_MUL_2X2(t3, t2, t1, t0, apre[1], apre[0], b[1], b[0]);
            FLINT_MPN_MUL_3P2X2(ucy, t1, t0, t3, t2, d[1], d[0]);
            FLINT_MPN_MUL_3P2X2(tcy, t3, t2, a[1], a[0], b[1], b[0]);
            sub_dddmmmsss(cy, t3, t2, tcy, t3, t2, ucy, t1, t0);
            /* We did an exact mulhigh, so one adjustment suffices. */
            if (cy || (t3 > d[1] || (t3 == d[1] && t2 >= d[0])))
                sub_ddmmss(t3, t2, t3, t2, d[1], d[0]);
            res[0] = t2;
            res[1] = t3;
        }
    }
    else
    {
        mp_ptr t, u;
        ulong tcy, ucy, cy;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(sizeof(mp_limb_t) * (2 * n));
        u = t + n;
        flint_mpn_mulhigh_n(t, apre, b, n);

        if (norm != 0)
        {
            flint_mpn_mullow_n(u, t, d, n);
            flint_mpn_mullow_n(t, a, b, n);
            mpn_sub_n(res, t, u, n);
            if (mpn_cmp(res, d, n) >= 0)
                mpn_sub_n(res, res, d, n);
        }
        else
        {
            /* Compute mod beta^(n+1) */
            ucy = flint_mpn_mullow_n(u, t, d, n);
            tcy = flint_mpn_mullow_n(t, a, b, n);
            cy = tcy - ucy - mpn_sub_n(res, t, u, n);
            /* We didn't do an exact mulhigh, so may need two adjustments */
            while (cy || mpn_cmp(res, d, n) >= 0)
                cy -= mpn_sub_n(res, res, d, n);
        }


        TMP_END;
    }
}

