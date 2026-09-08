/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* below this divisor length GMP's mpn_divisible_p is used when available */
#ifndef FLINT_MPN_DIVISIBLE_GMP_CUTOFF
#define FLINT_MPN_DIVISIBLE_GMP_CUTOFF 2048
#endif

/* from this dividend length on, a trial division of the divisor and the
   dividend by a primorial is tried first */
#ifndef FLINT_MPN_DIVISIBLE_PRIMORIAL_CUTOFF
#define FLINT_MPN_DIVISIBLE_PRIMORIAL_CUTOFF 32
#endif

#if FLINT_BITS != 64
/* 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23, for the 32-bit trial division */
#define PRIMORIAL UWORD(111546435)
#endif

/*
    Returns 1 if b divides a and 0 otherwise. Requires bn >= 1 and
    b[bn-1] != 0; a may have zero top limbs and an may be zero.

    Steps: the trivial cases (a = 0, |a| < |b|); 1 x 1 and 2 x 1 by
    hardware division; all other short inputs by GMP's mpn_divisible_p
    when available (its basecase has less overhead than the general path
    below); the 2-adic part (b = 2^v B^k b' with b' odd must divide the
    corresponding part of a); single-limb divisors via the Hensel remainder
    mpn_modexact_1_odd; for long dividends a cheap O(an) rejection by trial
    division: the residues
    of b and a modulo a primorial reveal a small prime dividing b but not
    a for roughly 70% of random pairs; then GMP's mpn_divisible_p for short
    divisors when available, and otherwise the Hensel division with
    remainder (flint_mpn_bdiv_qr) with n = an - bn + 1 quotient limbs:
    since 0 <= a, q b' < B^(n+bn), b' divides a iff the Hensel remainder
    (a - q b') / B^n mod B^bn vanishes.
*/
int
_flint_mpn_divisible(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t k = 0, n;
    unsigned int v;
    mp_ptr as, bs, q, r;
    int res;
    TMP_INIT;

    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    /* the smallest shapes by hardware division */
    if (bn == 1)
    {
        if (an == 1)
            return a[0] % b[0] == 0;
        if (an == 2)
        {
            if (a[1] == 0)
                return a[0] % b[0] == 0;
            mp_limb_t q, r0, r1;
            if (a[1] >= b[0])
            {
                r1 = a[1] % b[0];
                if (r1 == 0)
                    return a[0] % b[0] == 0;
            }
            else
                r1 = a[1];
            udiv_qrnnd(q, r0, r1, a[0], b[0]);
            (void) q;
            return r0 == 0;
        }
    }

#if FLINT_HAVE_NATIVE_mpn_divisible_p
    /* GMP's basecase has less overhead than anything below for short
       inputs; the trial-division filter only pays off for long dividends */
    if (an < FLINT_MPN_DIVISIBLE_PRIMORIAL_CUTOFF)
        return mpn_divisible_p(a, an, b, bn);
#endif

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an == 0)
        return 1;
    if (an < bn)
        return 0;

    /* zero low limbs of b */
    while (b[k] == 0)
        k++;
    if (k > 0)
    {
        if (!flint_mpn_zero_p(a, k))
            return 0;
        a += k;
        an -= k;
        b += k;
        bn -= k;
    }

    /* 2-adic valuation of b */
    v = flint_ctz(b[0]);
    if (v != 0 && (a[0] & ((UWORD(1) << v) - 1)) != 0)
        return 0;

    if (bn == 1)
    {
        mp_limb_t b0 = b[0] >> v;
        if (b0 == 1)
            return 1;
        /* b0 odd: b0 | a iff b0 | a / 2^v */
        return flint_mpn_divisible_1_odd(a, an, b0);
    }

    if (bn >= 4)
    {
#if FLINT_BITS == 64
        /* residues modulo 2^48 - 1 = 3^2 5 7 13 17 97 241 257 673 at 0.25
           ns per limb reveal a small prime dividing b but not a for about
           60% of random pairs */
        static const unsigned short ps[] = {3, 5, 7, 13, 17, 97, 241, 257, 673, 0};
        mp_limb_t rb = flint_mpn_mod_2exp48m1(b, bn) % ((UWORD(1) << 48) - 1), ra = 0;
#else
        static const unsigned short ps[] = {3, 5, 7, 11, 13, 17, 19, 23, 0};
        mp_limb_t rb = mpn_mod_1(b, bn, PRIMORIAL), ra = 0;
#endif
        int have_ra = 0, i;

        for (i = 0; ps[i] != 0; i++)
        {
            mp_limb_t p = ps[i];
            if (rb % p == 0)
            {
                if (!have_ra)
                {
#if FLINT_BITS == 64
                    ra = flint_mpn_mod_2exp48m1(a, an) % ((UWORD(1) << 48) - 1);
#else
                    ra = mpn_mod_1(a, an, PRIMORIAL);
#endif
                    have_ra = 1;
                }
                if (ra % p != 0)
                    return 0;
            }
        }
    }

#if FLINT_HAVE_NATIVE_mpn_divisible_p
    if (bn < FLINT_MPN_DIVISIBLE_GMP_CUTOFF)
        return mpn_divisible_p(a, an, b, bn);
#endif

    n = an - bn + 1;

    TMP_START;
    as = TMP_ALLOC((an + bn + n + bn) * sizeof(mp_limb_t));
    bs = as + an;
    q = bs + bn;
    r = q + n;

    if (v != 0)
    {
        mpn_rshift(as, a, an, v);
        mpn_rshift(bs, b, bn, v);
    }
    else
    {
        flint_mpn_copyi(as, a, an);
        flint_mpn_copyi(bs, b, bn);
    }

    flint_mpn_bdiv_qr(q, r, as, an, bs, bn, n);
    res = flint_mpn_zero_p(r, bn);

    TMP_END;
    return res;
}
