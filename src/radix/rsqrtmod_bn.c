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
    Reciprocal square root modulo B^n: develops y with

        x * y^2 == 1   (mod B^n),

    by Newton iteration, doubling the number of correct limbs each step exactly
    as radix_invmod_bn does for the reciprocal. The Newton update is

        y' = y (3 - x y^2) / 2 = y - y (x y^2 - 1) / 2.

    Writing x y^2 = 1 + B^m H (the low m limbs are 00..01 once y is correct to m
    limbs), this is y'[m, n) = -((y H)/2)[0, n-m).

    For ODD p, 2 is a unit, so (.)/2 mod B^k is the modular half (one high-to-low
    short division by 2, adding the modulus first when the value is odd), and the
    append-high-limbs structure above is exact: each step contributes the next
    n-m limbs and is never revisited.

    For p = 2 this structure fails: there 1/2 does not exist, the halving is an
    exact bit shift (y H >> 1, via radix_rshift_digits by one p-adic digit) that
    drops the top bit of the block it produces, and that lost bit sits at a low
    position which later steps never recompute. Instead, exactly as in
    _padic_sqrt_2, p = 2 runs a full-recompute Newton in bit precision: each step
    recomputes all of y' = y - y(x y^2 - 1)/2 to a new precision following the
    schedule e_{i+1} = (e_i + 3)/2 (each squaring gains 2k-2 bits), so the single
    bit lost by the shift is always at the advancing frontier.

    For odd p the high part H = (x y^2)[m, n) is obtained from a guarded middle
    product (its low m limbs are the identity 0..01), via the shared
    _radix_mulhigh_known_low, exactly as radix_invmod_bn forms (x y)[m, n). The
    p = 2 path still uses plain low products.

    Returns 1 on success, 0 if x is not a square modulo B (in which case nothing
    is written). y receives n limbs (caller provides room); y may alias x only
    if x has at least n limbs available (we read x while writing y), so callers
    should pass distinct buffers.
*/

/* keep only the low d bits of the limb array a (limbs hold e bits each, p = 2) */
static void
_radix2_mask(nn_ptr a, slong d, slong alimbs, slong e)
{
    slong full = d / e, r = d % e, i;

    if (r)
    {
        a[full] &= (UWORD(1) << r) - 1;
        for (i = full + 1; i < alimbs; i++)
            a[i] = 0;
    }
    else
    {
        for (i = full; i < alimbs; i++)
            a[i] = 0;
    }
}

/* W[0,k) = V[0,k) / 2 mod B^k, for ODD radix (2 invertible). One short division
   by two from the top, having (conceptually) added the modulus B^k when V is
   odd so that the dividend is even. W may alias V. */
static void
_radix_halve_odd(nn_ptr W, nn_srcptr V, slong k, const radix_t radix)
{
    ulong B = LIMB_RADIX(radix);
    ulong Bhalf = B >> 1;            /* (B - 1)/2, since B is odd */
    ulong rem, par;
    slong i;

    /* parity of V: B is odd, so V == sum of limbs (mod 2) */
    par = 0;
    for (i = 0; i < k; i++)
        par ^= V[i];
    par &= 1;

    /* short division of (V + par*B^k) by 2, from the top; avoids overflow that
       a literal rem*B + V[i] would cause when B > 2^63. */
    rem = par;
    for (i = k - 1; i >= 0; i--)
    {
        ulong v = V[i];
        if (rem)
        {
            W[i] = Bhalf + (v >> 1) + (v & 1);   /* (B + v) >> 1, no overflow */
            rem = 1 - (v & 1);                   /* (B + v) & 1, B odd */
        }
        else
        {
            W[i] = v >> 1;
            rem = v & 1;
        }
    }
    /* rem == 0 here: the dividend was even by construction */
}

/* single-limb (mod B = p^e) reciprocal square root; returns 0 if x is not a
   square unit modulo B. */
static ulong
n_rsqrtmod_bn(ulong x, const radix_t radix)
{
    nmod_t mod = radix->B;
    ulong p = DIGIT_RADIX(radix);
    slong e = radix->exp;
    ulong y, t;

    if (p == 2)
    {
        slong iter;

        if ((x & 7) != 1)
            return 0;                /* odd squares are 1 (mod 8) */

        y = 1;                       /* rsqrt of x mod 2^3 */
        /* Newton until x y^2 == 1 (mod B); converges in O(log e) steps, and the
           fixed point is stable, so a generous cap is safe. */
        for (iter = 0; iter < 2 * e; iter++)
        {
            t = nmod_mul(nmod_mul(y, y, mod), x, mod);   /* x y^2 mod B */
            if (t == 1)
                break;
            t = nmod_sub(t, 1, mod);                     /* even */
            t >>= 1;                                     /* (x y^2 - 1)/2 */
            y = nmod_sub(y, nmod_mul(y, t, mod), mod);
        }
        return y;
    }
    else
    {
        ulong inv2, s;
        slong d;

        s = n_sqrtmod(nmod_set_ui(x, radix->b), p);
        if (s == 0)
            return 0;                /* x is a unit but a non-residue mod p */

        /* canonical choice of residue */
        /* todo: it would be nice if this was simply enforced in n_sqrtmod */
        if (s > (p - 1) / 2)
            s = p - s;

        y = n_invmod(s, p);          /* rsqrt mod p */
        if (e == 1)
            return y;

        inv2 = (mod.n + 1) / 2;      /* 2^{-1} mod B (B odd) */
        for (d = 1; d < e; d *= 2)
        {
            t = nmod_mul(nmod_mul(y, y, mod), x, mod);   /* x y^2 mod B */
            t = nmod_sub(t, 1, mod);
            t = nmod_mul(t, inv2, mod);                  /* (x y^2 - 1)/2 */
            y = nmod_sub(y, nmod_mul(y, t, mod), mod);
        }
        return y;
    }
}

/* p = 2: develops y with x y^2 == 1 (mod B^n) by a full-recompute Newton in bit
   precision. res[0] must already hold the single-limb base (rsqrt mod 2^e). */
static void
_radix2_rsqrtmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)
{
    slong e = radix->exp;                 /* bits per limb (63 for p = 2) */
    slong D = e * n;                      /* total target bits */
    slong ed[2 * FLINT_BITS];
    slong L, k;
    ulong one = 1;
    nn_ptr y2, w, t;
    TMP_INIT;

    FLINT_ASSERT(e >= 3);

    /* zero-extend the base (correct mod 2^e) to n limbs */
    for (k = 1; k < n; k++)
        res[k] = 0;

    /* bit-precision schedule from D down to the base level (<= e) */
    ed[0] = D;
    L = 0;
    while (ed[L] > e)
    {
        if (L > 100) flint_abort();
        ed[L + 1] = (ed[L] + 3) / 2;
        L++;
    }

    TMP_START;
    y2 = TMP_ALLOC(n * sizeof(ulong));
    w  = TMP_ALLOC(n * sizeof(ulong));
    t  = TMP_ALLOC(n * sizeof(ulong));

    for (k = L - 1; k >= 0; k--)
    {
        slong nd = ed[k];                 /* target bits this step */
        slong hd = nd + 1;                /* one guard bit kept before the shift */
        slong nl = (nd + e - 1) / e;
        slong hl = (hd + e - 1) / e;

        /* w = (x y^2 - 1) mod 2^hd  (divisible by 2) */
        radix_mulmid(y2, res, hl, res, hl, 0, hl, radix);
        _radix2_mask(y2, hd, hl, e);
        radix_mulmid(w, x, FLINT_MIN(xn, hl), y2, hl, 0, hl, radix);
        _radix2_mask(w, hd, hl, e);
        radix_sub(w, w, hl, &one, 1, radix);
        _radix2_mask(w, hd, hl, e);

        /* w /= 2  ->  (x y^2 - 1)/2 mod 2^nd */
        radix_rshift_digits(w, w, hl, 1, radix);

        /* res = y - y w mod 2^nd */
        radix_mulmid(t, res, nl, w, nl, 0, nl, radix);
        _radix2_mask(t, nd, nl, e);
        radix_sub(res, res, nl, t, nl, radix);
        _radix2_mask(res, nd, nl, e);
    }

    TMP_END;
}

int
radix_rsqrtmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)
{
    slong i, m;
    slong a[FLINT_BITS];
    nn_ptr t, hbuf, scr;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= 1);

    res[0] = n_rsqrtmod_bn(x[0], radix);
    if (res[0] == 0)
        return 0;

    if (n == 1)
        return 1;

    if (DIGIT_RADIX(radix) == 2)
    {
        _radix2_rsqrtmod_bn(res, x, xn, n, radix);
        return 1;
    }

    /* odd p: Newton with the exact append-high-limbs structure */
    a[i = 0] = n;
    while (a[i] > 1)
    {
        a[i + 1] = (a[i] + 1) / 2;
        i++;
    }

    TMP_START;
    t   = TMP_ALLOC(n * sizeof(ulong));
    hbuf = TMP_ALLOC(n * sizeof(ulong));
    scr = TMP_ALLOC(n * sizeof(ulong));

    for (i--; i >= 0; i--)
    {
        slong nn = a[i];                 /* new size */
        slong nm;
        ulong one = 1;
        m = a[i + 1];                    /* current size (res correct to m limbs) */
        nm = nn - m;

        /* t = y^2 mod B^nn (square of the m-limb root). */
        radix_mulmid(t, res, m, res, m, 0, nn, radix);

        /* H = (x y^2)[m, nn): the low m limbs of x y^2 are the identity 00..01,
           so this is a high product with known low (a middle product with guard
           correction when m is large enough). */
        _radix_mulhigh_known_low(hbuf, x, FLINT_MIN(xn, nn), t, nn,
            &one, 1, m, nn, scr, radix);

        /* y'[m, nn) = -((y H)/2)[0, nm) */
        radix_mulmid(t, res, m, hbuf, nm, 0, nm, radix);     /* y H, low nm */
        _radix_halve_odd(t, t, nm, radix);                   /* (y H) 2^{-1} */
        radix_neg(res + m, t, nm, radix);
    }

    TMP_END;
    return 1;
}
