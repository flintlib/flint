/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

static void
radix_sqrmid_classical(nn_ptr z, nn_srcptr a, slong an, slong zlo, slong zhi, const radix_t radix)
{
    slong i, j;
    ulong hi, lo;
    ulong sbits;
    slong start, stop;
    nmod_t mod = radix->B;

    FLINT_ASSERT(an >= 1);

    an = FLINT_MIN(an, zhi);

    FLINT_ASSERT(zhi <= 2 * an);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zlo >= 0);

    sbits = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(an + 1);

    if (sbits <= FLINT_BITS)
    {
        ulong cy = 0;
        ulong s0 = 0;

        FLINT_ASSERT(mod.norm != 0);

        for (i = zlo; i < FLINT_MIN(zhi, 2 * an - 1); i++)
        {
            start = FLINT_MAX(0, i - an + 1);
            stop = FLINT_MIN(an - 1, (i + 1) / 2 - 1);

            s0 = 0;

            for (j = 0; j < stop - start + 1; j++)
            {
                s0 += a[start + j] * a[(i - stop) + (stop - start + 1) - 1 - j];
            }

            s0 = (s0 << 1);

            if (i % 2 == 0 && i / 2 < an)
                s0 += a[i / 2] * a[i / 2];

            cy += s0;

            z[i - zlo] = n_divrem_preinv_unnorm(&cy, cy, mod.n, mod.ninv, mod.norm);
        }

        if (zhi == 2 * an)
        {
            FLINT_ASSERT(cy < mod.n);
            z[2 * an - 1 - zlo] = cy;
        }
    }
    else if (sbits <= 2 * FLINT_BITS)
    {
        ulong cy[2] = { 0, 0, };
        ulong s0, s1;

        if (mod.norm == 0)
        {
            for (i = zlo; i < FLINT_MIN(zhi, 2 * an - 1); i++)
            {
                start = FLINT_MAX(0, i - an + 1);
                stop = FLINT_MIN(an - 1, (i + 1) / 2 - 1);

                s0 = s1 = 0;

                for (j = 0; j < stop - start + 1; j++)
                {
                    umul_ppmm(hi, lo, a[start + j], a[(i - stop) + (stop - start + 1) - 1 - j]);
                    add_ssaaaa(s1, s0, s1, s0, hi, lo);
                }

                s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
                s0 = (s0 << 1);

                if (i % 2 == 0 && i / 2 < an)
                {
                    umul_ppmm(hi, lo, a[i / 2], a[i / 2]);
                    add_ssaaaa(s1, s0, s1, s0, hi, lo);
                }

                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], s1, s0);

                z[i - zlo] = flint_mpn_divrem_2_1_preinv_norm(cy, cy, mod.n, mod.ninv);
            }
        }
        else
        {
            for (i = zlo; i < FLINT_MIN(zhi, 2 * an - 1); i++)
            {
                start = FLINT_MAX(0, i - an + 1);
                stop = FLINT_MIN(an - 1, (i + 1) / 2 - 1);

                s0 = s1 = 0;

                for (j = 0; j < stop - start + 1; j++)
                {
                    umul_ppmm(hi, lo, a[start + j], a[(i - stop) + (stop - start + 1) - 1 - j]);
                    add_ssaaaa(s1, s0, s1, s0, hi, lo);
                }

                s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
                s0 = (s0 << 1);

                if (i % 2 == 0 && i / 2 < an)
                {
                    umul_ppmm(hi, lo, a[i / 2], a[i / 2]);
                    add_ssaaaa(s1, s0, s1, s0, hi, lo);
                }

                add_ssaaaa(cy[1], cy[0], cy[1], cy[0], s1, s0);

                z[i - zlo] = flint_mpn_divrem_2_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }
        }

        if (zhi == 2 * an)
        {
            FLINT_ASSERT(cy[0] < mod.n);
            FLINT_ASSERT(cy[1] == 0);
            z[2 * an - 1 - zlo] = cy[0];
        }
    }
    else
    {
        if (mod.norm == 0)
        {
            ulong r;
            ulong cy0 = 0;
            ulong cy1 = 0;
            ulong cy2 = 0;
            ulong s0, s1, s2;

            for (i = zlo; i < FLINT_MIN(zhi, 2 * an - 1); i++)
            {
                start = FLINT_MAX(0, i - an + 1);
                stop = FLINT_MIN(an - 1, (i + 1) / 2 - 1);

                s0 = s1 = s2 = 0;

                for (j = 0; j + 3 < stop - start + 1; j += 4)
                {
                    umul_ppmm(hi, lo, a[start + j], a[(i - stop) + (stop - start + 1) - 1 - j]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[start + j + 1], a[(i - stop) + (stop - start + 1) - 1 - j - 1]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[start + j + 2], a[(i - stop) + (stop - start + 1) - 1 - j - 2]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[start + j + 3], a[(i - stop) + (stop - start + 1) - 1 - j - 3]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                }

                for ( ; j < stop - start + 1; j++)
                {
                    umul_ppmm(hi, lo, a[start + j], a[(i - stop) + (stop - start + 1) - 1 - j]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                }

                s2 = (s2 << 1) | (s1 >> (FLINT_BITS - 1));
                s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
                s0 = (s0 << 1);

                if (i % 2 == 0 && i / 2 < an)
                {
                    umul_ppmm(hi, lo, a[i / 2], a[i / 2]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                }

                add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, s2, s1, s0);

                /* cy2 is certainly already reduced */
                udiv_qrnnd_preinv(cy1, r, cy2, cy1, mod.n, mod.ninv);
                udiv_qrnnd_preinv(cy0, r, r, cy0, mod.n, mod.ninv);
                cy2 = 0;
                z[i - zlo] = r;
            }

            if (zhi == 2 * an)
            {
                FLINT_ASSERT(cy0 < mod.n);
                FLINT_ASSERT(cy1 == 0);
                FLINT_ASSERT(cy2 == 0);
                z[2 * an - 1 - zlo] = cy0;
            }
        }
        else
        {
            ulong cy[3] = { 0, 0, 0 };
            ulong s0, s1, s2;

            for (i = zlo; i < FLINT_MIN(zhi, 2 * an - 1); i++)
            {
                start = FLINT_MAX(0, i - an + 1);
                stop = FLINT_MIN(an - 1, (i + 1) / 2 - 1);

                s0 = s1 = s2 = 0;

                for (j = 0; j < stop - start + 1; j++)
                {
                    umul_ppmm(hi, lo, a[start + j], a[(i - stop) + (stop - start + 1) - 1 - j]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                }

                s2 = (s2 << 1) | (s1 >> (FLINT_BITS - 1));
                s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
                s0 = (s0 << 1);

                if (i % 2 == 0 && i / 2 < an)
                {
                    umul_ppmm(hi, lo, a[i / 2], a[i / 2]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
                }

                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], s2, s1, s0);

                z[i - zlo] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }

            if (zhi == 2 * an)
            {
                FLINT_ASSERT(cy[0] < mod.n);
                FLINT_ASSERT(cy[1] == 0);
                FLINT_ASSERT(cy[2] == 0);
                z[2 * an - 1 - zlo] = cy[0];
            }
        }
    }
}

void
radix_mulmid_classical(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    slong i, j, n1, n2;
    ulong hi, lo;
    ulong sbits;
    nmod_t mod = radix->B;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);

    if (an < bn)
    {
        FLINT_SWAP(nn_srcptr, a, b);
        FLINT_SWAP(slong, an, bn);
    }

    if (zhi <= 2 && zlo == 0)
    {
        if (zhi == 1)
        {
            z[0] = nmod_mul(a[0], b[0], radix->B);
            return;
        }
        else
        {
            ulong t[2];
            umul_ppmm(t[1], t[0], a[0], b[0]);
            z[0] = flint_mpn_divrem_1_preinv(t, t, 2, radix->B.n, radix->B.ninv, radix->B.norm);
            z[1] = t[0];
            if (bn >= 2)
            {
                z[1] = nmod_add(z[1], nmod_fmma(a[0], b[1], a[1], b[0], radix->B), radix->B);
            }
            else if (an >= 2)
            {
                z[1] = nmod_addmul(z[1], a[1], b[0], radix->B);
            }
            return;
        }
    }

    if (a == b && an == bn)
    {
        radix_sqrmid_classical(z, a, an, zlo, zhi, radix);
        return;
    }

    FLINT_ASSERT(an >= bn);

    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    FLINT_ASSERT(zhi <= an + bn);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zlo >= 0);

    /* Each dot product is bounded by bn*(B-1)^2 */
    /* Carry-in from lower terms is bounded by bn*B */
    /* Dot product with carry-in is bounded by (bn+1)*B^2 */

    sbits = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(bn + 1);

    if (sbits <= FLINT_BITS)
    {
        ulong cy = 0;

        FLINT_ASSERT(mod.norm != 0);

        for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
        {
            n1 = FLINT_MIN(an - 1, i);
            n2 = FLINT_MIN(bn - 1, i);

            for (j = 0; j < n1 + n2 - i + 1; j++)
                cy += a[i - n2 + j] * b[n2 - j];

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

        if (mod.norm == 0)
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                n1 = FLINT_MIN(an - 1, i);
                n2 = FLINT_MIN(bn - 1, i);

                for (j = 0; j < n1 + n2 - i + 1; j++)
                {
                    umul_ppmm(hi, lo, a[i - n2 + j], b[n2 - j]);
                    add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
                }

                z[i - zlo] = flint_mpn_divrem_2_1_preinv_norm(cy, cy, mod.n, mod.ninv);
            }
        }
        else
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                n1 = FLINT_MIN(an - 1, i);
                n2 = FLINT_MIN(bn - 1, i);

                for (j = 0; j < n1 + n2 - i + 1; j++)
                {
                    umul_ppmm(hi, lo, a[i - n2 + j], b[n2 - j]);
                    add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
                }

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
        if (mod.norm == 0)
        {
            ulong r;
            ulong cy0 = 0;
            ulong cy1 = 0;
            ulong cy2 = 0;

            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                n1 = FLINT_MIN(an - 1, i);
                n2 = FLINT_MIN(bn - 1, i);

                for (j = 0; j + 3 < n1 + n2 - i + 1; j += 4)
                {
                    umul_ppmm(hi, lo, a[i - n2 + j], b[n2 - j]);
                    add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[i - n2 + j + 1], b[n2 - j - 1]);
                    add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[i - n2 + j + 2], b[n2 - j - 2]);
                    add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                    umul_ppmm(hi, lo, a[i - n2 + j + 3], b[n2 - j - 3]);
                    add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                }

                for ( ; j < n1 + n2 - i + 1; j++)
                {
                    umul_ppmm(hi, lo, a[i - n2 + j], b[n2 - j]);
                    add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hi, lo);
                }

                /* cy2 is certainly already reduced */
                udiv_qrnnd_preinv(cy1, r, cy2, cy1, mod.n, mod.ninv);
                udiv_qrnnd_preinv(cy0, r, r, cy0, mod.n, mod.ninv);
                cy2 = 0;
                z[i - zlo] = r;
            }

            if (zhi == an + bn)
            {
                FLINT_ASSERT(cy0 < mod.n);
                FLINT_ASSERT(cy1 == 0);
                FLINT_ASSERT(cy2 == 0);
                z[an + bn - 1 - zlo] = cy0;
            }
        }
        else
        {
            ulong cy[3] = { 0, 0, 0 };

            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                n1 = FLINT_MIN(an - 1, i);
                n2 = FLINT_MIN(bn - 1, i);

                for (j = 0; j < n1 + n2 - i + 1; j++)
                {
                    umul_ppmm(hi, lo, a[i - n2 + j], b[n2 - j]);
                    add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], 0, hi, lo);
                }

                z[i - zlo] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);
            }

            if (zhi == an + bn)
            {
                FLINT_ASSERT(cy[0] < mod.n);
                FLINT_ASSERT(cy[1] == 0);
                FLINT_ASSERT(cy[2] == 0);
                z[an + bn - 1 - zlo] = cy[0];
            }
        }
    }
}

