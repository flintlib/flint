/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

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
    else if (sbits <= 2 * FLINT_BITS)
    {
        ulong cy[2] = { 0, 0, };

        if (mod.norm == 0)
        {
            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
                z[i] = flint_mpn_divrem_2_1_preinv_norm(cy, cy, mod.n, mod.ninv);
            }
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
                z[i] = flint_mpn_divrem_2_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }
        }

        FLINT_ASSERT(cy[0] < mod.n);
        FLINT_ASSERT(cy[1] == 0);
        return cy[0];
    }
    else
    {
        if (mod.norm == 0)
        {
            ulong r;
            ulong cy0 = 0;
            ulong cy1 = 0;
            ulong cy2 = 0;

            umul_ppmm(hi, lo, a[0], c);
            udiv_qrnnd_preinv(cy0, r, hi, lo, mod.n, mod.ninv);
            z[0] = r;

            for (i = 1; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                /* cy2 is certainly already reduced */
                udiv_qrnnd_preinv(cy1, r, cy2, cy1, mod.n, mod.ninv);
                udiv_qrnnd_preinv(cy0, r, r, cy0, mod.n, mod.ninv);
                cy2 = 0;
                z[i] = r;
            }

            FLINT_ASSERT(cy0 < mod.n);
            FLINT_ASSERT(cy1 == 0);
            return cy0;
        }
        else
        {
            ulong cy[3] = { 0, 0, 0 };

            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], 0, hi, lo);
                z[i] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }

            FLINT_ASSERT(cy[0] < mod.n);
            FLINT_ASSERT(cy[1] == 0);
            FLINT_ASSERT(cy[2] == 0);
            return cy[0];
        }
    }
}

