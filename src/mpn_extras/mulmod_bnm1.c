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
    r = a * b mod (B^rn - 1), as a representative < B^rn, following the
    splitting used by GMP's mulmod_bnm1: for even rn = 2h,

        u = a b mod (B^h - 1)      (recursively)
        v = a b mod (B^h + 1)      (flint_mpn_mulmod_2expp1_basecase)

    and, writing the product as p0 + p1 B^h, the reconstruction solves
    r ≡ u (mod B^h - 1), r ≡ v (mod B^h + 1) as

        r = u + (B^h - 1) s,   s = (u - v)/2 mod (B^h + 1),

    where the division by two adds B^h + 1 to odd values. Everything is
    a representative rather than canonical, which the recursion and the
    callers tolerate.
*/

/* even sizes at or below this use one plain multiplication: further
   splitting loses to the linear-pass overheads (residues, CRT, the
   2expp1 fold), which is also why the recursion stops at any odd
   size and why the rounding below never exceeds multiples of 8 */
#define MULMOD_BNM1_MUL_THRESHOLD 32

mp_size_t
flint_mpn_mulmod_bnm1_next_size(mp_size_t n)
{
    if (n <= MULMOD_BNM1_MUL_THRESHOLD)
        return n;
    if (n <= 96)
        return (n + 1) & ~(mp_size_t) 1;
    if (n <= 192)
        return (n + 3) & ~(mp_size_t) 3;
    return (n + 7) & ~(mp_size_t) 7;
}

mp_size_t
flint_mpn_mulmod_bnm1_itch(mp_size_t rn)
{
    return 8 * rn + 64;
}

/* rp (rn limbs) = value of (vp, vn) mod B^rn - 1, representative */
static void
bnm1_fold(nn_ptr rp, nn_srcptr vp, mp_size_t vn, mp_size_t rn)
{
    mp_size_t c = FLINT_MIN(vn, rn);
    ulong cy;

    flint_mpn_copyi(rp, vp, c);
    flint_mpn_zero(rp + c, rn - c);
    vp += c;
    vn -= c;

    while (vn > 0)
    {
        c = FLINT_MIN(vn, rn);
        cy = mpn_add(rp, rp, rn, vp, c);
        while (cy)
            cy = mpn_add_1(rp, rp, rn, cy);
        vp += c;
        vn -= c;
    }
}

void
flint_mpn_mulmod_bnm1(nn_ptr rp, mp_size_t rn, nn_srcptr ap, mp_size_t an,
                      nn_srcptr bp, mp_size_t bn, nn_ptr tp)
{
    mp_size_t h, ah, bh;
    nn_ptr a1, b1, a2, b2, v, s;
    int ca, cb, vtop, stop;
    ulong cy, br;

    FLINT_ASSERT(an >= 1 && bn >= 1);
    FLINT_ASSERT(an <= rn && bn <= rn);

    if (rn <= MULMOD_BNM1_MUL_THRESHOLD || (rn & 1))
    {
        if (ap == bp && an == bn)
            flint_mpn_sqr(tp, ap, an);
        else if (an >= bn)
            mpn_mul(tp, ap, an, bp, bn);
        else
            mpn_mul(tp, bp, bn, ap, an);
        bnm1_fold(rp, tp, an + bn, rn);
        return;
    }

    h = rn / 2;
    a1 = tp;
    b1 = a1 + h;
    a2 = b1 + h;
    b2 = a2 + h;
    v = b2 + h;          /* h limbs */
    s = v + h;           /* h + 1 limbs */
    tp = s + h + 1;      /* deeper levels; 2expp1 needs 2h here too */

    /* residues mod B^h - 1 (representative < B^h) and mod B^h + 1
       (value in [0, B^h]: h limbs plus a top flag) */
#define RESIDUES(x1, x2, cx, xp, xn) \
    do { \
        mp_size_t lo = FLINT_MIN(xn, h), hi = xn - lo; \
        cx = 0; \
        bnm1_fold(x1, xp, xn, h); \
        flint_mpn_copyi(x2, xp, lo); \
        flint_mpn_zero(x2 + lo, h - lo); \
        if (hi > 0) \
        { \
            br = mpn_sub(x2, x2, h, xp + h, hi); \
            if (br) \
            { \
                cy = mpn_add_1(x2, x2, h, 1);   /* += B^h + 1 */ \
                cx = (cy != 0);                 /* value B^h */ \
            } \
        } \
    } while (0)

    ah = an;
    bh = bn;
    RESIDUES(a1, a2, ca, ap, ah);
    if (ap == bp && an == bn)
    {
        b1 = a1;
        b2 = a2;
        cb = ca;
    }
    else
        RESIDUES(b1, b2, cb, bp, bh);
#undef RESIDUES

    /* u = a1 * b1 mod B^h - 1, into the low half of rp */
    flint_mpn_mulmod_bnm1(rp, h, a1, h, b1, h, tp);

    /* v = a2 * b2 mod B^h + 1, value vtop*B^h + v in [0, B^h] */
    vtop = flint_mpn_mulmod_2expp1_basecase(v, a2, b2,
                    ca | (cb << 1), FLINT_BITS * h, tp);

    /* s = (u - v)/2 mod B^h + 1, in [0, B^h]: h + 1 limbs */
    flint_mpn_copyi(s, rp, h);
    s[h] = 0;
    br = mpn_sub(s, s, h + 1, v, h);
    if (vtop)
        br += mpn_sub_1(s + h, s + h, 1, 1);
    if (br)
    {
        /* += B^h + 1 */
        cy = mpn_add_1(s, s, h + 1, 1);
        s[h] += 1;
        (void) cy;
    }
    if (s[0] & 1)
    {
        cy = mpn_add_1(s, s, h + 1, 1);
        s[h] += 1;
        (void) cy;
    }
    mpn_rshift(s, s, h + 1, 1);
    stop = (int) s[h];
    FLINT_ASSERT(stop <= 1);

    /* r = u + s * B^h - s  (mod B^rn - 1), u already in rp[0, h) */
    flint_mpn_copyi(rp + h, s, h);
    if (stop)                        /* s_h * B^rn ≡ s_h */
    {
        cy = mpn_add_1(rp, rp, rn, 1);
        while (cy)
            cy = mpn_add_1(rp, rp, rn, cy);
    }
    br = mpn_sub(rp, rp, rn, s, h + 1);
    while (br)                       /* -B^rn ≡ -1 */
        br = mpn_sub_1(rp, rp, rn, 1);
}
