/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "mpn_extras.h"
#include "fixed.h"

/*
    Integer square root with remainder by Newton-Karp-Markstein iteration,
    ported from radix_sqrtrem_newton_karp_markstein with B = 2^64.

    With sn = ceil(an/2), a is viewed as a fixed-point number alpha in
    [B^-2, 1) with 2 sn fraction limbs (zero-padded on top when an is odd),
    so that sqrt(a) = sqrt(alpha) B^sn. fixed_sqrt_newton with n2 = sn + 3
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
    fixed_sqrt_newton(S, Aview, viewn, n2);

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
    Two-limb square root using the hardware double-precision square root
    for the initial approximation (GMP's dedicated two-limb code does not).
    A double gives s to 53 bits: for a < 2^100 that is within one unit and
    only the final adjustment is needed, above that one Newton step
    s <- (s + a/s)/2 with a 128/64-bit division makes it exact up to one
    unit. Returns the number of remainder limbs (0, 1 or 2), writing the
    remainder to (r, 2) if r != NULL.
*/
static mp_size_t
_flint_mpn_sqrtrem_2(mp_ptr sp, mp_ptr r, mp_limb_t a1, mp_limb_t a0)
{
    mp_limb_t s, h, l, r1, r0;
    double d;

    d = (double) a1 * 18446744073709551616.0 + (double) a0;
    d = sqrt(d);
    if (d >= 18446744073709551615.0)
        s = UWORD_MAX;
    else
        s = (mp_limb_t) d;

    if (a1 >= (UWORD(1) << 36) && a1 != UWORD_MAX)
    {
        /* Newton step with a 128/64-bit division: a1 < s holds since
           s ~ sqrt(a) > a / 2^64 (enforced for the double's rounding) */
        mp_limb_t q, rem;
        if (s <= a1)
            s = a1 + 1;
        udiv_qrnnd(q, rem, a1, a0, s);
        s = s + (mp_limb_t) (((slong) (q - s)) / 2);
    }

    /* adjust so that s^2 <= a < (s+1)^2, tracking r = a - s^2 */
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
        /* stop when r < 2s + 1 */
        mp_limb_t t1, t0;
        t1 = s >> (FLINT_BITS - 1);
        t0 = (s << 1) + 1;
        t1 += (t0 == 0);
        if (r1 < t1 || (r1 == t1 && r0 < t0))
            break;
        sub_ddmmss(r1, r0, r1, r0, t1, t0);
        s++;
    }

    sp[0] = s;
    if (r != NULL)
    {
        r[0] = r0;
        r[1] = r1;
    }
    return (r1 != 0) ? 2 : (r0 != 0);
}
#endif

/*
    s = floor(sqrt(a)) with sn = ceil(an/2) limbs. If r != NULL, it is set
    to a - s^2 <= 2 s, zero-padded to sn + 1 limbs, and the number of limbs
    of the remainder is returned; r must have room for max(an, 2) limbs
    (GMP's mpn_sqrtrem, used below the Newton cutoff, needs an limbs of
    remainder space and writes in place). If r == NULL, the remainder is not formed
    and the return value is 0 for a perfect square and 1 otherwise (GMP's
    convention). Requires an >= 1 and a[an-1] != 0.
*/
mp_size_t
_flint_mpn_sqrtrem(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
{
    mp_size_t sn, rn;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    if (an == 1)
    {
        mp_limb_t r0;
        s[0] = n_sqrtrem(&r0, a[0]);
        if (r != NULL)
        {
            r[0] = r0;
            r[1] = 0;
        }
        return (r0 != 0);
    }

#if FLINT_BITS == 64
    if (an == 2)
    {
        rn = _flint_mpn_sqrtrem_2(s, r, a[1], a[0]);
        return (r == NULL) ? (rn != 0) : rn;
    }
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
