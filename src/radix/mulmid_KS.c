/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

void
radix_mulmid_KS(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    slong i;
    ulong sbits;
    nn_ptr tz, ta, tb;
    nmod_t mod = radix->B;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);

    if (an < bn)
    {
        FLINT_SWAP(nn_srcptr, a, b);
        FLINT_SWAP(slong, an, bn);
    }

    FLINT_ASSERT(an >= bn);

    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    FLINT_ASSERT(zhi <= an + bn);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zlo >= 0);

    TMP_START;

    sbits = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(bn + 1);

    if (sbits <= FLINT_BITS)
    {
        FLINT_ASSERT(mod.norm != 0);

        ulong cy = 0;

        ta = TMP_ALLOC((an + bn + (an + bn)) * sizeof(ulong));
        tb = ta + an;
        tz = tb + bn;

        for (i = 0; i < an; i++) ta[i] = a[i];
        for (i = 0; i < bn; i++) tb[i] = b[i];
        flint_mpn_mul(tz, ta, an, tb, bn);

        for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
        {
            cy += tz[i];
            z[i - zlo] = n_divrem_preinv_unnorm(&cy, cy, mod.n, mod.ninv, mod.norm);
        }

        if (zhi == an + bn)
        {
            FLINT_ASSERT(cy < mod.n);
            z[an + bn - 1 - zlo] = cy;
        }
    }
    else if (sbits <= 2 * FLINT_BITS)
    {
        ulong cy[2] = { 0, 0, };

        ta = TMP_ALLOC(2 * (an + bn + (an + bn)) * sizeof(ulong));
        tb = ta + 2 * an;
        tz = tb + 2 * bn;

        for (i = 0; i < an; i++) { ta[2 * i] = a[i]; ta[2 * i + 1] = 0; }
        for (i = 0; i < bn; i++) { tb[2 * i] = b[i]; tb[2 * i + 1] = 0; }
        flint_mpn_mul(tz, ta, 2 * an, tb, 2 * bn);

        if (mod.norm == 0)
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], tz[2 * i + 1], tz[2 * i]);
                z[i - zlo] = flint_mpn_divrem_2_1_preinv_norm(cy, cy, mod.n, mod.ninv);
            }
        }
        else
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], tz[2 * i + 1], tz[2 * i]);
                z[i - zlo] = flint_mpn_divrem_2_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }
        }

        if (zhi == an + bn)
        {
            FLINT_ASSERT(cy[0] < mod.n);
            FLINT_ASSERT(cy[1] == 0);
            z[an + bn - 1 - zlo] = cy[0];
        }
    }
    else
    {
        ulong cy[3] = { 0, 0, 0 };

        ta = TMP_ALLOC(3 * (an + bn + (an + bn)) * sizeof(ulong));
        tb = ta + 3 * an;
        tz = tb + 3 * bn;

        for (i = 0; i < an; i++) { ta[3 * i] = a[i]; ta[3 * i + 1] = 0; ta[3 * i + 2] = 0; }
        for (i = 0; i < bn; i++) { tb[3 * i] = b[i]; tb[3 * i + 1] = 0; tb[3 * i + 2] = 0; }
        flint_mpn_mul(tz, ta, 3 * an, tb, 3 * bn);

        if (mod.norm == 0)
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], tz[3 * i + 2], tz[3 * i + 1], tz[3 * i]);
                z[i - zlo] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, mod.n, mod.ninv);
            }
        }
        else
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], tz[3 * i + 2], tz[3 * i + 1], tz[3 * i]);
                z[i - zlo] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }
        }

        if (zhi == an + bn)
        {
            FLINT_ASSERT(cy[0] < mod.n);
            FLINT_ASSERT(cy[1] == 0);
            FLINT_ASSERT(cy[2] == 0);
            z[an + bn - 1 - zlo] = cy[0];
        }
    }

    TMP_END;
}

