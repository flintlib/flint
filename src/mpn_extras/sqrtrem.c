/*
    From GMP's mpn/generic/sqrtrem.c:
    Copyright 1999-2002, 2004, 2005, 2008, 2010, 2012, 2015, 2017 Free Software
    Foundation, Inc.
    Contributed to the GNU project by Paul Zimmermann (most code),
    Torbjorn Granlund (mpn_sqrtrem1) and Marco Bodrato (mpn_dc_sqrt).

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "mp_real.h"

/*
    Integer square root with remainder by Newton-Karp-Markstein iteration,
    ported from radix_sqrtrem_newton_karp_markstein with B = 2^64.

    With sn = ceil(an/2), a is viewed as a fixed-point number alpha in
    [B^-2, 1) with 2 sn fraction limbs (zero-padded on top when an is odd),
    so that sqrt(a) = sqrt(alpha) B^sn. _mp_real_sqrt_newton with n2 = sn + 3
    fraction limbs has absolute error at most 4 B^(-n2) / sqrt(alpha) <=
    4 B^(-n2+1), i.e. 4 B^-2 at the integer scale, so the integer part is
    certified when the first fraction limb lies in [2, B-2] and the
    remainder follows from a low product; otherwise the root is corrected
    by O(1) steps.

    s receives sn limbs. If r != NULL it receives a - s^2 <= 2 s,
    zero-padded to sn + 1 limbs. Requires an >= 2 and a[an-1] != 0.
*/
void
_flint_mpn_sqrtrem_newton(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
{
    mp_ptr S, P, t, q, Apad;
    mp_srcptr Aview;
    mp_size_t sn, n2, viewn;
    TMP_INIT;

    FLINT_ASSERT(an >= 2);
    FLINT_ASSERT(a[an - 1] != 0);

    sn = (an + 1) / 2;
    viewn = 2 * sn;

    TMP_START;

    if (an == viewn)
    {
        Aview = a;
    }
    else
    {
        Apad = TMP_ALLOC(viewn * sizeof(mp_limb_t));
        flint_mpn_copyi(Apad, a, an);
        Apad[viewn - 1] = 0;
        Aview = Apad;
    }

    n2 = sn + 3;

    S = TMP_ALLOC((n2 + 2) * sizeof(mp_limb_t));
    _mp_real_sqrt_newton(S, Aview, viewn, n2);

    FLINT_ASSERT(S[n2 + 1] == 0);
    /* q[0], ..., q[sn-1] are the candidate limbs of floor(sqrt(a)); q[-1]
       is the first fraction limb and q[sn] the integral overflow */
    q = S + n2 - sn;

    if (q[sn] == 0 && q[-1] > 1 && q[-1] < UWORD_MAX - 1)
    {
        if (r != NULL)
        {
            /* r = a - q^2 < 2 B^sn is determined by its low sn + 1 limbs
               (sn + 1 <= an since an >= 2) */
            P = TMP_ALLOC((sn + 1) * sizeof(mp_limb_t));
            flint_mpn_mulmid(P, q, sn, q, sn, 0, sn + 1);
            mpn_sub_n(r, a, P, sn + 1);
        }
        flint_mpn_copyi(s, q, sn);
        TMP_END;
        return;
    }

    if (q[sn] != 0)
        mpn_sub_1(q, q, sn + 1, 1);

    /* verification: P holds q^2, then a - q^2, over an + 1 limbs; t holds
       2q + 1 over sn + 1 limbs */
    P = TMP_ALLOC((an + 1 + sn + 1) * sizeof(mp_limb_t));
    t = P + an + 1;

    flint_mpn_mulmid(P, q, sn, q, sn, 0, FLINT_MIN(2 * sn, an + 1));
    if (2 * sn == an)
        P[an] = 0;      /* the only limb the product does not write */

    while (P[an] != 0 || mpn_cmp(P, a, an) > 0)
    {
        mpn_sub_1(q, q, sn, 1);
        /* q_old^2 - (2 q + 1) = q^2 */
        t[sn] = mpn_lshift(t, q, sn, 1);
        mpn_add_1(t, t, sn + 1, 1);
        mpn_sub(P, P, an + 1, t, sn + 1);
    }

    /* P = a - q^2 >= 0 */
    mpn_sub_n(P, a, P, an);
    P[an] = 0;

    for (;;)
    {
        t[sn] = mpn_lshift(t, q, sn, 1);
        mpn_add_1(t, t, sn + 1, 1);

        /* stop when a - q^2 < 2q + 1 */
        if (flint_mpn_zero_p(P + sn + 1, (an + 1) - (sn + 1))
            && mpn_cmp(P, t, sn + 1) < 0)
            break;

        mpn_add_1(q, q, sn, 1);
        mpn_sub(P, P, an + 1, t, sn + 1);
    }

    if (r != NULL)
        flint_mpn_copyi(r, P, sn + 1);
    flint_mpn_copyi(s, q, sn);

    TMP_END;
}

/* GMP fallback; r has room for an >= sn + 1 limbs, so mpn_sqrtrem can
   write its remainder (and use the space as scratch) in place */
void
_flint_mpn_sqrtrem_gmp(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
{
    if (r == NULL)
    {
        mpn_sqrtrem(s, NULL, a, an);
    }
    else
    {
        mp_size_t sn = (an + 1) / 2, rn;
        rn = mpn_sqrtrem(s, r, a, an);
        FLINT_ASSERT(rn <= sn + 1);
        flint_mpn_zero(r + rn, sn + 1 - rn);
    }
}

#if FLINT_BITS == 64

/*
    Two-limb square root. A double-precision square root gives the root to
    about 52 bits, which for a < 2^108 is within a few units, so that only
    the final adjustment is needed. Above that, one Newton step
    y = D + (a - D^2) / (2D) is taken from D ~ sqrt(a).

    With a fused multiply-add this step is evaluated entirely in double
    precision, which avoids the 128/64-bit division (slow on many x86 cores
    and done in software on ARM64). Working in units of 2^64, let
    d = sqrt(a1) rounded, so that D = 2^32 d is an integer when a1 >= 2^40,
    with |D - sqrt(a)| < 2^13; a0 only enters through the correction.
    Writing a / 2^64 = H + lo with H = a1 - (a1 mod 2^11) exact and
    lo = (a1 mod 2^11) + a0 / 2^64 (error < 2^-41), the correction is

        y - D = 2^31 (H - d^2 + lo) / d = t1 k + lo k,

    where t1 = H - d^2 (|t1| < 2^13) is computed by an FMA with a single
    rounding, and k = d q = 2^31 d / a1 ~ 2^31 / d (relative error ~2^-51)
    with q = 2^31 / a1 computed in parallel with the square root. The
    Newton step overestimates by less than 2^-27 and the rounding errors
    are below 2^-28 (2^-29 is observed), so after adding 2^-24 the
    truncation gives floor(sqrt(a)) <= s <= floor(sqrt(a)) + 1, where the
    upper value only occurs when sqrt(a) lies within ~2^-24 below an
    integer (the bias makes perfect squares come out exact; when d rounds
    up to 2^32, s may instead be up to 2 too small). The final adjustment
    makes the result exact in any case.

    Without a fast FMA, the Newton step uses a 128/64-bit division.
*/

/* fma() is a single instruction (clang defines no FP_FAST_FMA) */
#if defined(FP_FAST_FMA) || defined(__FP_FAST_FMA) || defined(__FMA__) \
    || defined(__ARM_FEATURE_FMA)
# define FLINT_SQRTREM_2_FMA 1
#endif

/* for a1 below this, the truncated double root is used as is */
#define FLINT_SQRTREM_2_SMALL (UWORD(1) << 44)

/* approximates floor(sqrt(a)) to a few units; requires a1 >= 2^40 */
static inline mp_limb_t
_sqrtrem_2_newton(mp_limb_t a1, mp_limb_t a0)
{
#ifdef FLINT_SQRTREM_2_FMA
    double a1d, a0d, q, d, H, lo, k, t1, t2, c;
    mp_limb_t s0;

    a1d = (double) a1;
    a0d = (double) a0;
    q = 0x1p31 / a1d;
    d = sqrt(a1d);
    H = (double) (a1 & ~UWORD(0x7ff));
    lo = fma(a0d, 0x1p-64, (double) (slong) (a1 & UWORD(0x7ff)));
    s0 = (d >= 0x1p32) ? UWORD_MAX : (mp_limb_t) (d * 0x1p32);
    k = d * q;
    t1 = fma(-d, d, H);
    /* the offset 16384 keeps c positive (|y - D| < 2^13) */
    t2 = fma(lo * q, d, 16384.0 + 0x1p-24);
    c = fma(t1, k, t2);
    return s0 + (mp_limb_t) (slong) c - 16384;
#else
    mp_limb_t s, q, rem;
    double d;

    d = sqrt((double) a1 * 0x1p64 + (double) a0);
    s = (d >= 0x1p64) ? UWORD_MAX : (mp_limb_t) d;

    /* for a1 = B - 1 the root is B - 1 */
    if (a1 != UWORD_MAX)
    {
        /* Newton step with a 128/64-bit division: a1 < s holds since
           s ~ sqrt(a) > a / 2^64 (enforced for the double's rounding) */
        if (s <= a1)
            s = a1 + 1;
        udiv_qrnnd(q, rem, a1, a0, s);
        s = s + (mp_limb_t) (((slong) (q - s)) / 2);
    }
    return s;
#endif
}

/* adjusts the estimate s so that s^2 <= a < (s+1)^2, setting
   (r1, r0) = a - s^2 */
static inline mp_limb_t
_sqrtrem_2_adjust(mp_limb_t * r1p, mp_limb_t * r0p, mp_limb_t s,
    mp_limb_t a1, mp_limb_t a0)
{
    mp_limb_t h, l, r1, r0, u1, u0;

    umul_ppmm(h, l, s, s);
    sub_ddmmss(r1, r0, a1, a0, h, l);
    while ((slong) r1 < 0)
    {
        /* s too large: r += 2s - 1, s -= 1 */
        s--;
        add_ssaaaa(r1, r0, r1, r0, 0, s);
        add_ssaaaa(r1, r0, r1, r0, 0, s);
        add_ssaaaa(r1, r0, r1, r0, 0, 1);
    }
    for (;;)
    {
        /* stop when r < 2s + 1; a single test on the sign of the
           difference, as comparing limbwise would branch unpredictably
           on r1 (which is 0 or 1 for large s) */
        sub_ddmmss(u1, u0, r1, r0, s >> (FLINT_BITS - 1), (s << 1) + 1);
        if ((slong) u1 < 0)
            break;
        r1 = u1;
        r0 = u0;
        s++;
    }
    *r1p = r1;
    *r0p = r0;
    return s;
}

/*
    The functions _flint_mpn_sqrtrem_n below follow the conventions of
    _flint_mpn_sqrtrem for an = n: if r != NULL, the remainder is written
    to (r, sn + 1) and its number of limbs returned; otherwise the return
    value is 0 for a perfect square and 1 otherwise.
*/
FLINT_STATIC_NOINLINE mp_size_t
_flint_mpn_sqrtrem_2(mp_ptr sp, mp_ptr r, mp_srcptr a)
{
    mp_limb_t a1 = a[1], a0 = a[0], s, r1, r0;

    if (a1 < FLINT_SQRTREM_2_SMALL)
#ifdef FLINT_SQRTREM_2_FMA
        s = (mp_limb_t) sqrt(fma((double) a1, 0x1p64, (double) a0));  /* exact product */
#else
        s = (mp_limb_t) sqrt((double) a1 * 0x1p64 + (double) a0);
#endif
    else
        s = _sqrtrem_2_newton(a1, a0);

    s = _sqrtrem_2_adjust(&r1, &r0, s, a1, a0);

    sp[0] = s;
    if (r == NULL)
        return (r1 | r0) != 0;
    r[0] = r0;
    r[1] = r1;
    return (r1 != 0) + ((r1 | r0) != 0);
}

/*
    Three- and four-limb square roots by one step of the divide and conquer
    algorithm of Zimmermann (Karatsuba square root) as in GMP's
    mpn_dc_sqrtrem, specialized to a two-limb root: for a normalized x
    (top limb >= 2^62) written as x = A beta^2 + x1 beta + x0 with
    x1, x0 < beta,

        (s', r') = sqrtrem(A)                      A has two limbs
        (q, u)   = divrem(r' beta + x1, 2 s')
        s = s' beta + q,  r = u beta + x0 - q^2
        if r < 0:  r += 2s - 1,  s -= 1

    For four limbs beta = B. For three limbs GMP appends a zero limb and
    discards the low 32 bits of the root; here the unbalanced split with
    beta = 2^32 is used instead (A is still the top two limbs), so that the
    division has a 33-bit quotient. The input is shifted left by 2c bits
    for normalization, and the root is shifted back as in GMP: with
    4^c a = S^2 + R and m = S mod 2^c, the root is S >> c and the
    remainder (R + m (2S - m)) / 4^c.
*/

/* s = floor(sqrt(a)) and (r1, r0) = a - s^2 <= 2s for a1 >= 2^62 */
static inline mp_limb_t
_sqrtrem_2_norm(mp_limb_t * r1, mp_limb_t * r0, mp_limb_t a1, mp_limb_t a0)
{
    return _sqrtrem_2_adjust(r1, r0, _sqrtrem_2_newton(a1, a0), a1, a0);
}

FLINT_STATIC_NOINLINE mp_size_t
_flint_mpn_sqrtrem_3(mp_ptr sp, mp_ptr rp, mp_srcptr a)
{
    mp_limb_t a2 = a[2], a1 = a[1], a0 = a[0];
    mp_limb_t x2, x1, x0, sh, rh1, rh0, n1, n0, Q, u, uc, q;
    mp_limb_t s1, s0, r1, r0, p1, p0, neg, m, w1, w0;
    unsigned int c, cc;

    /* branch-free normalization by 2c bits, 0 <= 2c <= 62 */
    c = flint_clz(a2) >> 1;
    cc = 2 * c;
    x2 = (a2 << cc) | ((a1 >> 1) >> (FLINT_BITS - 1 - cc));
    x1 = (a1 << cc) | ((a0 >> 1) >> (FLINT_BITS - 1 - cc));
    x0 = (a0 << cc);

    /* 2^63 <= sh < 2^64, rh <= 2 sh */
    sh = _sqrtrem_2_norm(&rh1, &rh0, x2, x1);

    /* N = rh 2^32 + (x0 >> 32), N / sh < 2^34 */
    n1 = (rh1 << 32) | (rh0 >> 32);
    n0 = (rh0 << 32) | (x0 >> 32);
    udiv_qrnnd(Q, u, n1, n0, sh);
    /* q = N div 2 sh <= 2^32, (uc, u) = N mod 2 sh */
    q = Q >> 1;
    add_ssaaaa(uc, u, 0, u, 0, sh & (-(Q & 1)));

    /* s = sh 2^32 + q */
    add_ssaaaa(s1, s0, sh >> 32, sh << 32, 0, q);

    /* r = u 2^32 + (x0 mod 2^32) - q^2, |r| < 2^98 */
    umul_ppmm(p1, p0, q, q);
    sub_ddmmss(r1, r0, (uc << 32) | (u >> 32), (u << 32) | (x0 & UWORD(0xffffffff)), p1, p0);

    /* if r < 0: s -= 1, r += 2s + 1 (without branching: this happens
       with probability ~1/4) */
    neg = -(r1 >> (FLINT_BITS - 1));
    sub_ddmmss(s1, s0, s1, s0, 0, neg & 1);
    add_ssaaaa(r1, r0, r1, r0, neg & ((s1 << 1) | (s0 >> (FLINT_BITS - 1))), neg & ((s0 << 1) | 1));

    /* undo the normalization */
    m = s0 & ((UWORD(1) << c) - 1);
    sp[0] = (s0 >> c) | ((s1 << 1) << (FLINT_BITS - 1 - c));
    sp[1] = s1 >> c;

    if (rp == NULL)
        return (r1 | r0 | m) != 0;

    /* 4^c r = R + m (2S - m) < 2^128, with 2S - m < 2^97 */
    sub_ddmmss(w1, w0, (s1 << 1) | (s0 >> (FLINT_BITS - 1)), s0 << 1, 0, m);
    umul_ppmm(p1, p0, m, w0);
    add_ssaaaa(r1, r0, r1, r0, p1 + m * w1, p0);

    r0 = (r0 >> cc) | ((r1 << 1) << (FLINT_BITS - 1 - cc));
    r1 = r1 >> cc;
    rp[0] = r0;
    rp[1] = r1;
    rp[2] = 0;

    /* r1 = 0 is rare here (r < 2^97 is usually large) */
    return (r1 != 0) ? 2 : (r0 != 0);
}

/* s = (s1, s0) = floor(sqrt(x)) and x - s^2 = (r2, r1, r0) <= 2s (r2 is
   returned) for a normalized x = (x3, x2, x1, x0), x3 >= 2^62 */
static inline mp_limb_t
_sqrtrem_4_norm(mp_limb_t * s1p, mp_limb_t * s0p, mp_limb_t * r1p, mp_limb_t * r0p,
    mp_limb_t x3, mp_limb_t x2, mp_limb_t x1, mp_limb_t x0)
{
    mp_limb_t sh, rh1, rh0, qh, ql, u, uc, t, s1, s0, r2, r1, r0, p1, p0, neg;

    sh = _sqrtrem_2_norm(&rh1, &rh0, x3, x2);

    /* divide N = rh B + x1 < (2 sh + 1) B by sh, as GMP: first reduce the
       top limb below sh (at most twice), counting in qh */
    qh = rh1;
    rh0 -= sh & (-rh1);
    t = (rh0 >= sh);
    rh0 -= sh & (-t);
    qh += t;
    udiv_qrnnd(ql, u, rh0, x1, sh);

    /* N div sh = qh B + ql = 2 s_low + (ql & 1); N mod 2 sh = u + (ql & 1) sh.
       s_low = (qh >> 1) B + s0 <= B, so s0 = 0 when qh >> 1 = 1 */
    s0 = (ql >> 1) | (qh << (FLINT_BITS - 1));
    qh >>= 1;
    add_ssaaaa(uc, u, 0, u, 0, sh & (-(ql & 1)));

    /* r = (uc, u) B + x0 - s_low^2 */
    umul_ppmm(p1, p0, s0, s0);
    sub_dddmmmsss(r2, r1, r0, uc, u, x0, qh, p1, p0);

    /* s = sh B + s_low; s1 wraps to 0 if sh = B - 1 and s_low = B, in
       which case r < 0 and the correction unwraps it */
    s1 = sh + qh;
    neg = -(r2 >> (FLINT_BITS - 1));
    sub_ddmmss(s1, s0, s1, s0, 0, neg & 1);
    add_sssaaaaaa(r2, r1, r0, r2, r1, r0, neg & (s1 >> (FLINT_BITS - 1)),
        neg & ((s1 << 1) | (s0 >> (FLINT_BITS - 1))), neg & ((s0 << 1) | 1));

    *s1p = s1;
    *s0p = s0;
    *r1p = r1;
    *r0p = r0;
    return r2;
}

FLINT_STATIC_NOINLINE mp_size_t
_flint_mpn_sqrtrem_4(mp_ptr sp, mp_ptr rp, mp_srcptr a)
{
    mp_limb_t a3 = a[3], a2 = a[2], a1 = a[1], a0 = a[0];
    mp_limb_t x3, x2, x1, x0, t, s1, s0, r2, r1, r0, p1, p0, m, w2, w1, w0;
    unsigned int c, cc;

    c = flint_clz(a3) >> 1;
    cc = 2 * c;
    x3 = (a3 << cc) | ((a2 >> 1) >> (FLINT_BITS - 1 - cc));
    x2 = (a2 << cc) | ((a1 >> 1) >> (FLINT_BITS - 1 - cc));
    x1 = (a1 << cc) | ((a0 >> 1) >> (FLINT_BITS - 1 - cc));
    x0 = (a0 << cc);

    r2 = _sqrtrem_4_norm(&s1, &s0, &r1, &r0, x3, x2, x1, x0);

    m = s0 & ((UWORD(1) << c) - 1);
    sp[0] = (s0 >> c) | ((s1 << 1) << (FLINT_BITS - 1 - c));
    sp[1] = s1 >> c;

    if (rp == NULL)
        return (r2 | r1 | r0 | m) != 0;

    /* 4^c r = R + m (2S - m) */
    sub_dddmmmsss(w2, w1, w0, s1 >> (FLINT_BITS - 1), (s1 << 1) | (s0 >> (FLINT_BITS - 1)), s0 << 1, 0, 0, m);
    umul_ppmm(p1, p0, m, w0);
    umul_ppmm(w1, t, m, w1);
    add_ssaaaa(w1, t, w1, t, 0, p1);
    add_sssaaaaaa(r2, r1, r0, r2, r1, r0, w1 + (m & (-w2)), t, p0);

    r0 = (r0 >> cc) | ((r1 << 1) << (FLINT_BITS - 1 - cc));
    r1 = (r1 >> cc) | ((r2 << 1) << (FLINT_BITS - 1 - cc));
    r2 = r2 >> cc;
    rp[0] = r0;
    rp[1] = r1;
    rp[2] = r2;

    return (r2 != 0) + ((r2 | r1) != 0) + ((r2 | r1 | r0) != 0);
}


/*
    Divide and conquer (Karatsuba) square root with remainder, as GMP's
    mpn_dc_sqrtrem: for {np, 2n} with np[2n-1] >= 2^62, sets {sp, n} to
    the square root and {np, n} + c B^n to the remainder (<= 2s), where the
    carry c is returned. With l = floor(n/2), h = n - l the root of the top
    2h limbs is computed recursively, the next l limbs of the root are the
    quotient of (r' B^l + next limbs) by 2 s', and the remainder is
    corrected once (the division is _flint_mpn_divrem_preinv1). The recursion
    stops at h <= 2, where the two- and four-limb code above is used. The limbs {np + n, n} are destroyed.
    scratch needs n + 2 limbs.

    If approx != 0 (only at the top level, where approx = 2^k - 2 when the
    low k bits of the root will be discarded), the function may return 2
    early with only the root computed: when (sp[0] & approx) != 0, the low
    k bits of the root are >= 2, so the final correction cannot change the
    root after discarding them, and the input cannot be a square.
*/
static mp_limb_t
_sqrtrem_divconquer_norm(mp_ptr sp, mp_ptr np, mp_size_t n, mp_limb_t approx, mp_ptr scratch)
{
    mp_limb_t q, b;
    slong c;
    mp_size_t l, h;

    /* the callers (with an >= 5 input limbs) and the recursion have n >= 2 */
    FLINT_ASSERT(n >= 2);

    if (n == 2)
        return _sqrtrem_4_norm(sp + 1, sp, np + 1, np, np[3], np[2], np[1], np[0]);

    l = n / 2;
    h = n - l;

    q = _sqrtrem_divconquer_norm(sp + l, np + 2 * l, h, 0, scratch);
    if (q != 0)
        mpn_sub_n(np + 2 * l, np + 2 * l, sp + l, h);

    /* {np + l, n} / {sp + l, h}: quotient in {scratch, l} plus the high
       limb, remainder in {np + l, h} */
    q += _flint_mpn_divrem_preinv1(scratch, np + l, n, sp + l, h,
        flint_mpn_preinv1(sp[n - 1], sp[n - 2]), scratch + l + 1);
    c = scratch[0] & 1;
    mpn_rshift(sp, scratch, l, 1);
    sp[l - 1] |= q << (FLINT_BITS - 1);
    if (FLINT_UNLIKELY((sp[0] & approx) != 0))
        return 2;
    q >>= 1;

    if (c != 0)
        c = mpn_add_n(np + l, np + l, sp + l, h);

    flint_mpn_sqr(np + n, sp, l);
    b = q + mpn_sub_n(np, np, np + n, 2 * l);
    c -= (l == h) ? (slong) b : (slong) mpn_sub_1(np + 2 * l, np + 2 * l, 1, b);

    if (c < 0)
    {
        q = mpn_add_1(sp + l, sp + l, h, q);
        c += mpn_addmul_1(np, sp, n, 2) + 2 * q;
        c -= mpn_sub_1(np, np, n, 1);
        q -= mpn_sub_1(sp, sp, n, 1);
    }

    return c;
}

/*
    Square root without remainder, a port of GMP's mpn_dc_sqrt (Marco
    Bodrato): writes in {sp, n} the square root (rounded towards zero) of
    {np, 2n - odd}, and returns 0 if the operand is a perfect square and
    1 otherwise. Requires {np, 2n - odd} 4^nsh to be normalized, n > 4
    and nsh < FLINT_BITS / 2.

    With l = (n - 1) / 2 < h = n - l, the top 2h limbs (shifted) are passed
    to _sqrtrem_divconquer_norm, and the next root limbs come from an
    approximate division with one extra limb, which (as h > l) also
    absorbs the q^2 term of the remainder; only when that extra limb is
    too small to be conclusive is the remainder evaluated.
*/
static int
_sqrt_divconquer(mp_ptr sp, mp_srcptr np, mp_size_t n, unsigned int nsh, unsigned int odd)
{
    mp_limb_t q, cy;
    int c;
    mp_size_t l, h;
    mp_ptr qp, tp, scratch;
    TMP_INIT;

    FLINT_ASSERT(np[2 * n - 1 - odd] != 0);
    FLINT_ASSERT(n > 4);
    FLINT_ASSERT(nsh < FLINT_BITS / 2);

    l = (n - 1) / 2;
    h = n - l;

    TMP_START;
    /* scratch: n + 1 limbs; tp: n + h + 2 limbs (the quotient
       {tp + n + 1, l + 2} overlaps its top), with tp[-1] writable */
    scratch = TMP_ALLOC((2 * n + h + 4) * sizeof(mp_limb_t));
    tp = scratch + n + 2;

    if (nsh != 0)
    {
        /* o is used to exactly set the lowest bits of the dividend */
        int o = l > (1 + odd);
        mpn_lshift(tp - o, np + l - 1 - o - odd, n + h + 1 + o, 2 * nsh);
    }
    else
    {
        flint_mpn_copyi(tp, np + l - 1 - odd, n + h + 1);
    }

    q = _sqrtrem_divconquer_norm(sp + l, tp + l + 1, h, 0, scratch);
    if (q != 0)
        mpn_sub_n(tp + l + 1, tp + l + 1, sp + l, h);

    qp = tp + n + 1;   /* l + 2 limbs */
    flint_mpn_divapprox(qp, tp, n + 1, sp + l, h);
    q += qp[l + 1];
    c = 1;

    if (q > 1)
    {
        flint_mpn_store(sp, l, UWORD_MAX);
    }
    else
    {
        mpn_rshift(sp, qp + 1, l, 1);
        sp[l - 1] |= q << (FLINT_BITS - 1);

        if (((qp[0] >> 3) | (qp[1] & (UWORD_MAX >> ((FLINT_BITS >> odd) - nsh - 1)))) == 0)
        {
            /* The approximation is not good enough: the extra limb (and
               nsh bits) is smaller than the possible error.
               {qp + 1, l + 1} equals 2 {sp, l}. The remainder
               R = {tp + 1, n} - {sp + l, h} {qp + 1, l + 1} of the
               approximate division satisfies -2 s' <= R < B^h (GMP computes
               the full product; unlike GMP, only its low h + 1 limbs are
               computed here, which determine R). */
            flint_mpn_mulmid(scratch, sp + l, h, qp + 1, l + 1, 0, h + 1);
            mpn_sub_n(tp + 1, tp + 1, scratch, h + 1);
            if ((slong) tp[1 + h] < 0)
            {
                /* only if the approximate quotient was too large */
                cy = mpn_addmul_1(tp + 1, sp + l, h, 2);
                tp[1 + h] += cy;
                mpn_sub_1(sp, sp, l, 1);
            }

            if (tp[1 + h] == 0 && flint_mpn_zero_p(tp + l + 1, h - l))
            {
                flint_mpn_sqr(scratch, sp, l);
                c = mpn_cmp(tp + 1, scratch + l, l);
                if (c == 0)
                {
                    if (nsh != 0)
                    {
                        mpn_lshift(tp, np, l, 2 * nsh);
                        np = tp;
                    }
                    c = mpn_cmp(np, scratch + odd, l - odd);
                }
                if (c < 0)
                {
                    mpn_sub_1(sp, sp, l, 1);
                    c = 1;
                }
            }
        }
    }

    TMP_END;

    if ((odd | nsh) != 0)
        mpn_rshift(sp, sp, n, nsh + (odd ? FLINT_BITS / 2 : 0));

    return c;
}

/*
    Square root by _sqrtrem_divconquer_norm, for an >= 5. As in GMP, the input
    is shifted left by 2c bits and, when an is odd, padded with a zero limb
    at the bottom, so that 4^k a = S^2 + R with k = c (+ FLINT_BITS / 2);
    then s = S >> k and r = (R + m (2S - m)) / 4^k with m = S mod 2^k.
*/
mp_size_t
_flint_mpn_sqrtrem_divconquer(mp_ptr sp, mp_ptr rp, mp_srcptr a, mp_size_t an)
{
    mp_size_t tn, sn, rn, i;
    mp_ptr tp, scratch;
    mp_limb_t rl, m, cy;
    unsigned int c, k;
    int odd;
    TMP_INIT;

    FLINT_ASSERT(an >= 5);
    FLINT_ASSERT(a[an - 1] != 0);

    sn = tn = (an + 1) / 2;
    odd = an & 1;
    c = flint_clz(a[an - 1]) / 2;
    k = c + (odd ? FLINT_BITS / 2 : 0);

    if (rp == NULL && an > 8)
        return _sqrt_divconquer(sp, a, tn, c, odd);

    TMP_START;
    tp = TMP_ALLOC((3 * tn + 3) * sizeof(mp_limb_t));
    scratch = tp + 2 * tn + 1;

    tp[0] = 0;
    if (c != 0)
        mpn_lshift(tp + odd, a, an, 2 * c);
    else
        flint_mpn_copyi(tp + odd, a, an);

    rl = _sqrtrem_divconquer_norm(sp, tp, tn, (rp == NULL && k != 0) ? (UWORD(1) << k) - 2 : 0, scratch);

    if (rl == 2)
    {
        /* early exit: not a square, root correct after the shift */
        mpn_rshift(sp, sp, tn, k);
        TMP_END;
        return 1;
    }

    if (k != 0)
    {
        /* R += 2 m S - m^2, S >>= k */
        m = sp[0] & ((UWORD(1) << k) - 1);
        rl += mpn_addmul_1(tp, sp, tn, 2 * m);
        cy = mpn_submul_1(tp, &m, 1, m);
        rl -= mpn_sub_1(tp + 1, tp + 1, tn - 1, cy);
        mpn_rshift(sp, sp, tn, k);
    }

    if (rp == NULL)
    {
        rn = (rl != 0) || !flint_mpn_zero_p(tp, tn);
        TMP_END;
        return rn;
    }

    /* r = {tp, tn + 1} >> 2k */
    tp[tn] = rl;
    k *= 2;
    rn = tn + 1;
    if (k >= FLINT_BITS)
    {
        tp++;
        rn--;
        k -= FLINT_BITS;
    }
    if (k != 0)
        mpn_rshift(rp, tp, rn, k);
    else
        flint_mpn_copyi(rp, tp, rn);
    for (i = rn; i < sn + 1; i++)
        rp[i] = 0;

    TMP_END;

    rn = sn + 1;
    while (rn > 0 && rp[rn - 1] == 0)
        rn--;
    return rn;
}
#endif

FLINT_STATIC_NOINLINE mp_size_t
_flint_mpn_sqrtrem_large(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an);

/*
    s = floor(sqrt(a)) with sn = ceil(an/2) limbs. If r != NULL, it is set
    to a - s^2 <= 2 s, zero-padded to sn + 1 limbs, and the number of limbs
    of the remainder is returned; r must have room for max(an, 2) limbs
    (GMP's mpn_sqrtrem, used below the Newton cutoff on 32-bit machines,
    needs an limbs of remainder space and writes in place). If r == NULL, the remainder is not formed
    and the return value is 0 for a perfect square and 1 otherwise (GMP's
    convention). Requires an >= 1 and a[an-1] != 0.
*/
mp_size_t
_flint_mpn_sqrtrem(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
{
    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    if (an == 1)
    {
        /* as n_sqrtrem */
        mp_limb_t s0, r0;
        s0 = (mp_limb_t) sqrt((double) a[0]);
        s0 -= (s0 * s0 > a[0]);
#if FLINT_BITS == 64
        if (FLINT_UNLIKELY(s0 == (UWORD(1) << 32)))
            s0--;
#endif
        r0 = a[0] - s0 * s0;
        s[0] = s0;
        if (r != NULL)
        {
            r[0] = r0;
            r[1] = 0;
        }
        return (r0 != 0);
    }

#if FLINT_BITS == 64
    /* separate functions, so that this dispatch needs no stack frame */
    if (an == 2)
        return _flint_mpn_sqrtrem_2(s, r, a);

    if (an == 3)
        return _flint_mpn_sqrtrem_3(s, r, a);

    if (an == 4)
        return _flint_mpn_sqrtrem_4(s, r, a);
#endif

    return _flint_mpn_sqrtrem_large(s, r, a, an);
}

/* kept out of line so that the small cases do not pay for its stack frame */
FLINT_STATIC_NOINLINE mp_size_t
_flint_mpn_sqrtrem_large(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
{
    mp_size_t sn, rn;

#if FLINT_BITS == 64
    if (an < FLINT_MPN_SQRTREM_NEWTON_CUTOFF)
        return _flint_mpn_sqrtrem_divconquer(s, r, a, an);
#endif

    sn = (an + 1) / 2;

    if (r == NULL)
    {
        /* only the root and the exactness are wanted */
        if (an < FLINT_MPN_SQRTREM_NEWTON_CUTOFF)
            return mpn_sqrtrem(s, NULL, a, an) != 0;

        _flint_mpn_sqrtrem_newton(s, NULL, a, an);

        /* a = s^2 is tested first on the low two limbs, then in full */
        {
            mp_limb_t h, l;
            umul_ppmm(h, l, s[0], s[0]);
            h += 2 * s[0] * s[1];
            if (l != a[0] || h != a[1])
                return 1;
        }

        {
            mp_ptr t;
            int exact;
            TMP_INIT;
            TMP_START;
            t = TMP_ALLOC(2 * sn * sizeof(mp_limb_t));
            flint_mpn_sqr(t, s, sn);
            exact = (mpn_cmp(t, a, an) == 0) && (an == 2 * sn || t[2 * sn - 1] == 0);
            TMP_END;
            return !exact;
        }
    }

    if (an < FLINT_MPN_SQRTREM_NEWTON_CUTOFF)
        _flint_mpn_sqrtrem_gmp(s, r, a, an);
    else
        _flint_mpn_sqrtrem_newton(s, r, a, an);

    rn = sn + 1;
    while (rn > 0 && r[rn - 1] == 0)
        rn--;

    return rn;
}
