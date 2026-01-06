/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

FLINT_FORCE_INLINE mp_limb_t
flint_mpn_divrem_3_1_preinv_norm_inline(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
{
    mp_limb_t r;
    r = n_divrem_norm(&qp[2], up[2], d);
    udiv_qrnnd_preinv(qp[1], r, r, up[1], d, dinv);
    udiv_qrnnd_preinv(qp[0], r, r, up[0], d, dinv);
    return r;
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
        ulong cy[3] = { 0, 0, 0 };

        if (mod.norm == 0)
        {
            for (i = zlo; i < FLINT_MIN(zhi, an + bn - 1); i++)
            {
                n1 = FLINT_MIN(an - 1, i);
                n2 = FLINT_MIN(bn - 1, i);

#define cy0 cy[0]
#define cy1 cy[1]
#define cy2 cy[2]

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

#undef cy0
#undef cy1
#undef cy2

                z[i - zlo] = flint_mpn_divrem_3_1_preinv_norm_inline(cy, cy, mod.n, mod.ninv);
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
                    add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], 0, hi, lo);
                }

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
}

