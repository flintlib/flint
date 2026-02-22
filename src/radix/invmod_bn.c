/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "radix.h"

/* TODO: consider storing a dynamic table in the radix context object. */

static const uint8_t inv2_tab[] = { 0, 1 };
static const uint8_t inv3_tab[] = { 0, 1, 2 };
static const uint8_t inv4_tab[] = { 0, 1, 0, 3 };
static const uint8_t inv5_tab[] = { 0, 1, 3, 2, 4 };
static const uint8_t inv6_tab[] = { 0, 1, 0, 0, 0, 5 };
static const uint8_t inv7_tab[] = { 0, 1, 4, 5, 2, 3, 6 };
static const uint8_t inv8_tab[] = { 0, 1, 0, 3, 0, 5, 0, 7 };
static const uint8_t inv9_tab[] = { 0, 1, 5, 0, 7, 2, 0, 4, 8 };
static const uint8_t inv10_tab[] = { 0, 1, 0, 7, 0, 0, 0, 3, 0, 9 };
static const uint8_t inv11_tab[] = { 0, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10 };
static const uint8_t inv12_tab[] = { 0, 1, 0, 0, 0, 5, 0, 7, 0, 0, 0, 11 };
static const uint8_t inv13_tab[] = { 0, 1, 7, 9, 10, 8, 11, 2, 5, 3, 4, 6, 12 };
static const uint8_t inv14_tab[] = { 0, 1, 0, 5, 0, 3, 0, 0, 0, 11, 0, 9, 0, 13 };
static const uint8_t inv15_tab[] = { 0, 1, 8, 0, 4, 0, 0, 13, 2, 0, 0, 11, 0, 7, 14 };
static const uint8_t inv16_tab[] = { 0, 1, 0, 11, 0, 13, 0, 7, 0, 9, 0, 3, 0, 5, 0, 15 };
static const uint8_t inv17_tab[] = { 0, 1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16 };
static const uint8_t inv18_tab[] = { 0, 1, 0, 0, 0, 11, 0, 13, 0, 0, 0, 5, 0, 7, 0, 0, 0, 17 };
static const uint8_t inv19_tab[] = { 0, 1, 10, 13, 5, 4, 16, 11, 12, 17, 2, 7, 8, 3, 15, 14, 6, 9, 18 };
static const uint8_t inv20_tab[] = { 0, 1, 0, 7, 0, 0, 0, 3, 0, 9, 0, 11, 0, 17, 0, 0, 0, 13, 0, 19 };
static const uint8_t inv21_tab[] = { 0, 1, 11, 0, 16, 17, 0, 0, 8, 0, 19, 2, 0, 13, 0, 0, 4, 5, 0, 10, 20 };
static const uint8_t inv22_tab[] = { 0, 1, 0, 15, 0, 9, 0, 19, 0, 5, 0, 0, 0, 17, 0, 3, 0, 13, 0, 7, 0, 21 };
static const uint8_t inv23_tab[] = { 0, 1, 12, 8, 6, 14, 4, 10, 3, 18, 7, 21, 2, 16, 5, 20, 13, 19, 9, 17, 15, 11, 22 };
static const uint8_t inv24_tab[] = { 0, 1, 0, 0, 0, 5, 0, 7, 0, 0, 0, 11, 0, 13, 0, 0, 0, 17, 0, 19, 0, 0, 0, 23 };
static const uint8_t inv25_tab[] = { 0, 1, 13, 17, 19, 0, 21, 18, 22, 14, 0, 16, 23, 2, 9, 0, 11, 3, 7, 4, 0, 6, 8, 12, 24 };
static const uint8_t inv26_tab[] = { 0, 1, 0, 9, 0, 21, 0, 15, 0, 3, 0, 19, 0, 0, 0, 7, 0, 23, 0, 11, 0, 5, 0, 17, 0, 25 };
static const uint8_t inv27_tab[] = { 0, 1, 14, 0, 7, 11, 0, 4, 17, 0, 19, 5, 0, 25, 2, 0, 22, 8, 0, 10, 23, 0, 16, 20, 0, 13, 26 };
static const uint8_t inv28_tab[] = { 0, 1, 0, 19, 0, 17, 0, 0, 0, 25, 0, 23, 0, 13, 0, 15, 0, 5, 0, 3, 0, 0, 0, 11, 0, 9, 0, 27 };
static const uint8_t inv29_tab[] = { 0, 1, 15, 10, 22, 6, 5, 25, 11, 13, 3, 8, 17, 9, 27, 2, 20, 12, 21, 26, 16, 18, 4, 24, 23, 7, 19, 14, 28 };

static ulong
n_invmod_bn(ulong x, unsigned int exp, nmod_t bmod, nmod_t mod)
{
    ulong r, y;
    unsigned int e;

    FLINT_ASSERT(bmod.n >= 2);

    switch (bmod.n)
    {
        case 2: r = inv2_tab[x % 2]; break;
        case 3: r = inv3_tab[x % 3]; break;
        case 4: r = inv4_tab[x % 4]; break;
        case 5: r = inv5_tab[x % 5]; break;
        case 6: r = inv6_tab[x % 6]; break;
        case 7: r = inv7_tab[x % 7]; break;
        case 8: r = inv8_tab[x % 8]; break;
        case 9: r = inv9_tab[x % 9]; break;
        case 10: r = inv10_tab[x % 10]; break;
        case 11: r = inv11_tab[x % 11]; break;
        case 12: r = inv12_tab[x % 12]; break;
        case 13: r = inv13_tab[x % 13]; break;
        case 14: r = inv14_tab[x % 14]; break;
        case 15: r = inv15_tab[x % 15]; break;
        case 16: r = inv16_tab[x % 16]; break;
        case 17: r = inv17_tab[x % 17]; break;
        case 18: r = inv18_tab[x % 18]; break;
        case 19: r = inv19_tab[x % 19]; break;
        case 20: r = inv20_tab[x % 20]; break;
        case 21: r = inv21_tab[x % 21]; break;
        case 22: r = inv22_tab[x % 22]; break;
        case 23: r = inv23_tab[x % 23]; break;
        case 24: r = inv24_tab[x % 24]; break;
        case 25: r = inv25_tab[x % 25]; break;
        case 26: r = inv26_tab[x % 26]; break;
        case 27: r = inv27_tab[x % 27]; break;
        case 28: r = inv28_tab[x % 28]; break;
        case 29: r = inv29_tab[x % 29]; break;
        default:
            y = n_gcdinv(&r, nmod_set_ui(x, bmod), mod.n);
            if (y != 1)
                r = 0;
    }

    if (exp >= 2 && r != 0)
    {
        y = nmod_mul(x, r, mod) - 1;
        r = nmod_mul(r, nmod_sub(1, y, mod), mod);

        for (e = 2; e < exp; e *= 2)
        {
            y = nmod_mul(y, y, mod);
            r = nmod_mul(r, y + 1, mod);
        }
    }

    return r;
}

int
radix_invmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)
{
    res[0] = n_invmod_bn(x[0], radix->exp, radix->b, radix->B);
    if (res[0] == 0)
        return 0;

    if (n == 1)
        return 1;

    ulong t2[2];
    ulong u2[2];
    ulong v2[2];

    /* 1 -> 2 */
    /* mul_1x1 */
    umul_ppmm(t2[1], t2[0], res[0], res[0]);
    u2[0] = flint_mpn_divrem_1_preinv(t2, t2, 2, radix->B.n, radix->B.ninv, radix->B.norm);
    u2[1] = t2[0];
    radix_mulmid(v2, u2, 2, x, FLINT_MIN(xn, 2), 0, 2, radix);
    res[1] = n_negmod(v2[1], radix->B.n);

    if (n == 2)
        return 1;

    nn_ptr u;
    TMP_INIT;
    TMP_START;
    /* Todo: can tighten allocation when we always do a mulhigh */
    u = TMP_ALLOC(n * sizeof(ulong));

    slong a[FLINT_BITS];
    slong i, m;
    a[i = 0] = n;
    while (n > 2)
        a[++i] = (n = (n + 1) / 2);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        slong rxn = FLINT_MIN(n, m + xn);

        /* Can use middle product */
        if (m > 3 && LIMB_RADIX(radix) >= m && rxn > (m - 3))
        {
            ulong one = 1;
            radix_mulmid(u + m - 3, res, m, x, FLINT_MIN(n, xn), m - 3, rxn, radix);
            /* We know that the m least significant limbs of the full product
               res * x are 00...001. The approximate high product is either
               exact or 1 ulp too small. The latter is indicated by a
               nonzero m-1 limb (i.e. a borrow to be returned). */
            if (u[m - 1] != 0)
                radix_add(u + m, u + m, n - m, &one, 1, radix);
            radix_mulmid(res + m, u + m, rxn - m, res, m, 0, n - m, radix);
        }
        else
        {
            radix_mulmid(u, res, m, x, FLINT_MIN(n, xn), 0, rxn, radix);
            radix_mulmid(res + m, u + m, rxn - m, res, m, 0, n - m, radix);
        }
        radix_neg(res + m, res + m, n - m, radix);
    }

    TMP_END;

    return 1;
}

