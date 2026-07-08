/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Windowed middle product via a balanced high product.  flint_mpn_mulhigh_n
    returns (a lower approximation of) the top n limbs, i.e. limbs [n, 2n), of a
    balanced n x n product.

    Only the high ends of a and b reach the window: a limb a[p] can contribute to
    limb >= zlo only if p + (bc - 1) >= zlo, and a limb at or above zhi never
    reaches the window at all.  So we first clamp to zhi and slice off the low,
    non-contributing limbs

        ac = min(an, zhi),  bc = min(bn, zhi),
        pa = max(0, zlo - (bc - 1)),  qb = max(0, zlo - (ac - 1)),
        a' = a + pa  (La = ac - pa limbs),  b' = b + qb  (Lb = bc - qb limbs),
        zlo2 = zlo - pa - qb,

    which leaves a balanced-ish problem of size ~ zn rather than ~ max(an, bn).
    The sliced product a'*b' still contains a few sub-zlo diagonals, but they lie
    below limb zlo2 and are dropped by the high product just like the ones sliced
    away -- so the net effect is the required "drop p + q < zlo".

    The window [zlo2, zlo2 + zn) of a'*b' is then exposed as a top slice: pad both
    operands to n limbs and, when the window starts below M = max(La, Lb), shift
    them up by L = M - zlo2 low zero limbs so that limb zlo2 lands at limb n:

        L   = max(0, M - zlo2),        n = M + L,        off = zlo2 - M + L (>= 0),
        res = mulhigh_n(0^L,a',0.. , 0^L,b',0..)  ~=  limbs [n,2n) of a'*b'*B^{2L},
        z   = res[off .. off + zn).

    Correct for any input; the result is a lower approximation of the exact
    window (mulhigh never exceeds the true high part).  Cheapest when the window
    sits at the top of the product (zlo2 >= M, L = 0).

    Copies are avoided where possible: when L = 0 an operand of length n is passed
    verbatim and only a shorter one is padded; padded operands zero only their
    genuine zero limbs; and when the window is exactly mulhigh_n's whole output
    (off = 0, zn = n) it is written straight into z.
*/
void
flint_mpn_mulmid_via_mulhigh_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                               mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t zn = zhi - zlo;
    mp_size_t ac, bc, pa, qb, La, Lb, zlo2, M, L, n, off;
    mp_srcptr xp, yp;
    mp_ptr res;
    TMP_INIT;

    FLINT_ASSERT(an >= 1 && bn >= 1);
    FLINT_ASSERT(0 <= zlo && zlo < zhi && zhi <= an + bn);

    ac = FLINT_MIN(an, zhi);
    bc = FLINT_MIN(bn, zhi);
    pa = (zlo > bc - 1) ? zlo - (bc - 1) : 0;
    qb = (zlo > ac - 1) ? zlo - (ac - 1) : 0;
    La = ac - pa;
    Lb = bc - qb;

    if (La <= 0 || Lb <= 0)
    {
        /* every partial product lands below zlo: the lower-approximation window
           is all zero */
        flint_mpn_zero(z, zn);
        return;
    }

    a += pa;
    b += qb;
    zlo2 = zlo - pa - qb;

    M = FLINT_MAX(La, Lb);
    L = (zlo2 < M) ? (M - zlo2) : 0;
    n = M + L;
    off = (zlo2 - M) + L;
    xp = a;
    yp = b;

    FLINT_ASSERT(off >= 0 && off + zn <= n);

    TMP_START;

    if (L == 0)
    {
        if (La < n)
        {
            mp_ptr X = TMP_ARRAY_ALLOC(n, mp_limb_t);
            flint_mpn_copyi(X, a, La);
            flint_mpn_zero(X + La, n - La);
            xp = X;
        }
        if (Lb < n)
        {
            mp_ptr Y = TMP_ARRAY_ALLOC(n, mp_limb_t);
            flint_mpn_copyi(Y, b, Lb);
            flint_mpn_zero(Y + Lb, n - Lb);
            yp = Y;
        }
    }
    else
    {
        mp_ptr X = TMP_ARRAY_ALLOC(n, mp_limb_t);
        mp_ptr Y = TMP_ARRAY_ALLOC(n, mp_limb_t);

        flint_mpn_zero(X, L);
        flint_mpn_copyi(X + L, a, La);
        flint_mpn_zero(X + L + La, n - L - La);

        flint_mpn_zero(Y, L);
        flint_mpn_copyi(Y + L, b, Lb);
        flint_mpn_zero(Y + L + Lb, n - L - Lb);

        xp = X;
        yp = Y;
    }

    if (off == 0 && zn == n)
    {
        flint_mpn_mulhigh_n(z, xp, yp, n);
    }
    else
    {
        res = TMP_ARRAY_ALLOC(n, mp_limb_t);
        flint_mpn_mulhigh_n(res, xp, yp, n);
        flint_mpn_copyi(z, res + off, zn);
    }

    TMP_END;
}
