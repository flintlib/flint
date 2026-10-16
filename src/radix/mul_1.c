/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* Inline version of flint_mpn_divrem_2_1_preinv_unnorm for speed.  */
FLINT_FORCE_INLINE mp_limb_t
flint_mpn_divrem_2_1_preinv_unnorm1(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)
{
    mp_limb_t u0, u1, r;

    FLINT_ASSERT(norm >= 1);

    u1 = up[1];
    u0 = up[0];
    if (u1 < d)
    {
        d <<= norm;
        qp[1] = 0;
        r = (u1 << norm) | (u0 >> (FLINT_BITS - norm));
    }
    else
    {
        /* Hack: this branch is unreachable, but leaving it in results
           in faster code with GCC on Zen 3. */
        d <<= norm;
        r = (u1 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    }

    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}


ulong
radix_mul_1(nn_ptr z, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    slong i;
    ulong hi, lo;
    slong sbits;
    nmod_t mod = radix->B;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(c < mod.n);

    sbits = 2 * NMOD_BITS(mod);

    if (sbits <= FLINT_BITS)
    {
        ulong cy = 0;
        FLINT_ASSERT(mod.norm != 0);

        for (i = 0; i < n; i++)
        {
            cy += a[i] * c;
            z[i] = n_divrem_preinv_unnorm(&cy, cy, mod.n, mod.ninv, mod.norm);
        }

        FLINT_ASSERT(cy < mod.n);
        return cy;
    }
    else
    {
        if (mod.norm == 0)
        {
            ulong cy0 = 0;
            ulong cy1 = 0;
            ulong r;

            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_ssaaaa(cy1, cy0, cy1, cy0, hi, lo);
                r = n_divrem_norm(&cy1, cy1, mod.n);
                udiv_qrnnd_preinv(cy0, r, r, cy0, mod.n, mod.ninv);
                z[i] = r;
            }

            FLINT_ASSERT(cy0 < mod.n);
            FLINT_ASSERT(cy1 == 0);
            return cy0;
        }
        else
        {
            ulong cy[2] = { 0, 0, };

            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
                z[i] = flint_mpn_divrem_2_1_preinv_unnorm1(cy, cy, mod.n, mod.ninv, mod.norm);
            }

            FLINT_ASSERT(cy[0] < mod.n);
            FLINT_ASSERT(cy[1] == 0);
            return cy[0];
        }
    }
}

