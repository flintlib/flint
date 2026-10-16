/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/*
    Fused add/sub with an on-the-fly left shift by sh p-adic digits, 1 <= sh < e.
    These fuse the carry loop of radix_add / radix_sub with the within-limb shift
    of radix_lshift_digits, so the shifted operand never has to be materialised in
    a temporary. The shift is restricted to a single limb's worth of digits
    (1 <= sh < e); the radix_integer wrappers decompose a general shift into a limb
    offset plus such a digit shift.

    The shifted operand (b << sh) spans bn + 1 limbs:
        limb 0      : (b[0] mod be2) * be
        limb i      : (b[i] mod be2) * be + (b[i-1] div be2)        (0 < i < bn)
        limb bn     : (b[bn-1] div be2)
    where be = p^sh and be2 = p^{e-sh}. Each shifted limb is a valid limb in
    [0, B): the low part (b[i] mod be2) * be < B and the carried high part
    (b[i-1] div be2) < be, and their sum is < B.

    In every routine res may alias a (a[i] is read before res[i] is written at the
    same index), but res must NOT alias b (b is read while res is written).
*/

/* hi = x / be2, *lo = x mod be2, via the radix's precomputed division. */
#define RADIX_LSH_DIVREM(hi, lo, x, be2, pre)                       \
    do {                                                            \
        if ((pre)->m == 0)                                          \
            (hi) = n_divrem_precomp_m0(&(lo), (x), (be2), (pre));   \
        else if ((pre)->c == 0)                                     \
            (hi) = n_divrem_precomp_c0(&(lo), (x), (be2), (pre));   \
        else                                                        \
            (hi) = n_divrem_precomp_c1_bounded(&(lo), (x), (be2), (pre)); \
    } while (0)

/* res = a + (b << sh).  Returns the carry out of limb max(an, bn+1) (0 or 1). */
ulong
radix_addlsh(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn,
    unsigned int sh, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    ulong be = radix->bpow[sh];
    ulong be2 = radix->bpow[radix->exp - sh];
    n_div_precomp_t pre;
    ulong shcy = 0, sb, qh, rl;
    slong i, L;

    FLINT_ASSERT(sh >= 1 && sh < radix->exp);
    FLINT_ASSERT(bn >= 1);

    *pre = radix->bpow_div[radix->exp - sh];
    L = FLINT_MAX(an, bn + 1);

    cy = 0;
    cy -= 1;            /* inverted carry, exactly as radix_add */

    for (i = 0; i < L; i++)
    {
        if (i < bn)
        {
            RADIX_LSH_DIVREM(qh, rl, b[i], be2, pre);
            sb = shcy + rl * be;
            shcy = qh;
        }
        else if (i == bn)
        {
            sb = shcy;
        }
        else
        {
            sb = 0;
        }

        {
            ulong ai = (i < an) ? a[i] : 0;
            sub_ddmmss(hi, lo, 0, ai + (cy + 1), 0, B - sb);
            res[i] = lo + (hi & B);
            cy = hi;
        }
    }

    cy += 1;
    return cy;
}

/* res = a - (b << sh).  Requires a >= (b << sh); returns the borrow (0). */
ulong
radix_sublsh(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn,
    unsigned int sh, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    ulong be = radix->bpow[sh];
    ulong be2 = radix->bpow[radix->exp - sh];
    n_div_precomp_t pre;
    ulong shcy = 0, sb, qh, rl;
    slong i, L;

    FLINT_ASSERT(sh >= 1 && sh < radix->exp);
    FLINT_ASSERT(bn >= 1);

    *pre = radix->bpow_div[radix->exp - sh];
    L = FLINT_MAX(an, bn + 1);

    cy = 0;

    for (i = 0; i < L; i++)
    {
        if (i < bn)
        {
            RADIX_LSH_DIVREM(qh, rl, b[i], be2, pre);
            sb = shcy + rl * be;
            shcy = qh;
        }
        else if (i == bn)
        {
            sb = shcy;
        }
        else
        {
            sb = 0;
        }

        {
            ulong ai = (i < an) ? a[i] : 0;
            sub_ddmmss(hi, lo, 0, ai, 0, sb - cy);
            res[i] = lo + (hi & B);
            cy = hi;
        }
    }

    return -cy;          /* 0 when a >= b << sh */
}

/* res = (b << sh) - a.  Requires (b << sh) >= a; returns the borrow (0). */
ulong
radix_lshsub(nn_ptr res, nn_srcptr b, slong bn, unsigned int sh,
    nn_srcptr a, slong an, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    ulong be = radix->bpow[sh];
    ulong be2 = radix->bpow[radix->exp - sh];
    n_div_precomp_t pre;
    ulong shcy = 0, sb, qh, rl;
    slong i, L;

    FLINT_ASSERT(sh >= 1 && sh < radix->exp);
    FLINT_ASSERT(bn >= 1);

    *pre = radix->bpow_div[radix->exp - sh];
    L = FLINT_MAX(an, bn + 1);

    cy = 0;

    for (i = 0; i < L; i++)
    {
        if (i < bn)
        {
            RADIX_LSH_DIVREM(qh, rl, b[i], be2, pre);
            sb = shcy + rl * be;
            shcy = qh;
        }
        else if (i == bn)
        {
            sb = shcy;
        }
        else
        {
            sb = 0;
        }

        {
            ulong ai = (i < an) ? a[i] : 0;
            sub_ddmmss(hi, lo, 0, sb, 0, ai - cy);
            res[i] = lo + (hi & B);
            cy = hi;
        }
    }

    return -cy;          /* 0 when b << sh >= a */
}

/* sign of |a| - |b << sh|; a and b assumed normalised (no leading zero limbs). */
int
radix_cmplsh(nn_srcptr a, slong an, nn_srcptr b, slong bn, unsigned int sh,
    const radix_t radix)
{
    ulong be = radix->bpow[sh];
    ulong be2 = radix->bpow[radix->exp - sh];
    n_div_precomp_t pre;
    ulong qh, rl, sblo, sbhi, sb, alimb;
    slong byn, j;

    FLINT_ASSERT(sh >= 1 && sh < radix->exp);

    if (bn == 0)
        return (an == 0) ? 0 : 1;
    if (an == 0)
        return -1;

    *pre = radix->bpow_div[radix->exp - sh];

    /* length of b << sh: bn, or bn+1 if the top sh digits of b[bn-1] are set */
    RADIX_LSH_DIVREM(qh, rl, b[bn - 1], be2, pre);
    byn = bn + (qh != 0);

    if (an != byn)
        return (an < byn) ? -1 : 1;

    /* equal length: compare from the top */
    for (j = byn - 1; j >= 0; j--)
    {
        /* low part of shifted limb j: (b[j] mod be2) * be */
        if (j < bn)
        {
            RADIX_LSH_DIVREM(qh, rl, b[j], be2, pre);
            sblo = rl * be;
        }
        else
            sblo = 0;

        /* high part carried in: (b[j-1] div be2) */
        if (j - 1 >= 0 && j - 1 < bn)
        {
            RADIX_LSH_DIVREM(sbhi, rl, b[j - 1], be2, pre);
        }
        else
            sbhi = 0;

        sb = sblo + sbhi;
        alimb = a[j];
        if (alimb != sb)
            return (alimb < sb) ? -1 : 1;
    }

    return 0;
}
