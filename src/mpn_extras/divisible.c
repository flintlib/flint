/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* a trial division of the divisor and the dividend by small primes is
   tried first when both the divisor and the quotient have at least this
   many limbs; for divisible inputs it costs a few percent there */
#define DIVISIBLE_SCREEN_CUTOFF 16

#if FLINT_BITS != 64
/* 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23, for the 32-bit trial division */
#define PRIMORIAL UWORD(111546435)
#endif

/*
    Returns 1 if b divides a and 0 otherwise. Requires bn >= 1 and
    b[bn-1] != 0; a may have zero top limbs and an may be zero.

    Steps: the trivial cases (a = 0, |a| < |b|); 1 x 1 and 2 x 1 by
    hardware division; the 2-adic part (b = 2^v B^k b' with b' odd must
    divide the corresponding part of a); divisors b' of one limb via the
    Hensel remainder mpn_modexact_1_odd; for divisors and quotients of at
    least DIVISIBLE_SCREEN_CUTOFF limbs a cheap O(an) rejection by trial
    division: the residues of b and a modulo small primes reveal a prime
    dividing b but not a for most random pairs; finally the Hensel
    remainder as in GMP (_flint_mpn_divisible_bdiv), or for huge operands
    the Newton-based exact division q = a / b mod B^n with n = an - bn + 1,
    which is the quotient if b divides a, and a comparison of q b with a.
*/
int
_flint_mpn_divisible(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t k = 0, n, qn;
    unsigned int v;
    mp_ptr q, t;
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

    if (bn == 2 && v != 0 && (b[1] >> v) == 0)
    {
        /* b / 2^v fits in one limb */
        mp_limb_t b0 = (b[0] >> v) | (b[1] << (FLINT_BITS - v));
        return flint_mpn_divisible_1_odd(a, an, b0);
    }

    /* short odd parts of the divisor: the Hensel remainder in registers */
    if (bn - ((b[bn - 1] >> v) == 0) <= FLINT_MPN_DIVEXACT_SMALL_BN)
        return _flint_mpn_divisible_small(a, an, b, bn, v);

    if (bn >= DIVISIBLE_SCREEN_CUTOFF && an - bn >= DIVISIBLE_SCREEN_CUTOFF)
    {
        mp_limb_t rb, ra, mask, amask;

        /* residues modulo 2^48 - 1 = 3^2 5 7 13 17 97 241 257 673 at 0.25
           ns per limb (respectively modulo 3 5 7 11 13 17 19 23 on 32-bit
           machines) reveal a small prime dividing b but not a for about
           60% of random pairs; mask gets the primes dividing b, amask those
           not dividing a (the moduli being constants, the reductions are
           multiplications) */
#if FLINT_BITS == 64
# define SCREEN_PRIMES(X) X(0, 3) X(1, 5) X(2, 7) X(3, 13) X(4, 17) X(5, 97) X(6, 241) X(7, 257) X(8, 673)
# define SCREEN_RES(x, xn) (flint_mpn_mod_2exp48m1(x, xn) % ((UWORD(1) << 48) - 1))
#else
# define SCREEN_PRIMES(X) X(0, 3) X(1, 5) X(2, 7) X(3, 11) X(4, 13) X(5, 17) X(6, 19) X(7, 23)
# define SCREEN_RES(x, xn) mpn_mod_1(x, xn, PRIMORIAL)
#endif
#define SCREEN_BMASK(i, p) mask |= (mp_limb_t) (rb % (p) == 0) << (i);
#define SCREEN_AMASK(i, p) amask |= (mp_limb_t) (ra % (p) != 0) << (i);

        rb = SCREEN_RES(b, bn);
        mask = 0;
        SCREEN_PRIMES(SCREEN_BMASK)

        if (mask != 0)
        {
            ra = SCREEN_RES(a, an);
            amask = 0;
            SCREEN_PRIMES(SCREEN_AMASK)
            if ((mask & amask) != 0)
                return 0;
        }

#undef SCREEN_PRIMES
#undef SCREEN_RES
#undef SCREEN_BMASK
#undef SCREEN_AMASK
    }

    n = an - bn + 1;

    /* the Hensel remainder as in GMP's mpn_divisible_p, except for huge
       operands where the Newton-based exact division is faster */
    if (!(FLINT_MIN(n, bn) >= FLINT_MPN_DIVEXACT_NEWTON_CUTOFF
            || (bn >= FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF && n >= 8 * bn)))
        return _flint_mpn_divisible_bdiv(a, an, b, bn, v);

    /* b divides a iff q b = a for the quotient q = a / b mod B^qn given by
       the exact division, where qn = n, or n - 1 when the top limb of a is
       below that of b (the exact division then sets q[n - 1] = 0). The
       low qn limbs of q b agree with a by construction, so only the high
       limbs of the product are formed (a middle product corrected with
       the known low limbs) and compared with the high limbs of a. */
    qn = n - (n > 1 && a[an - 1] < b[bn - 1]);

    TMP_START;
    q = TMP_ALLOC((n + (an + 1 - qn) + an + 1) * sizeof(mp_limb_t));
    t = q + n;

    _flint_mpn_divexact(q, a, an, b, bn);

    if (qn >= bn)
        _flint_mpn_mulhigh_known_low(t, q, qn, b, bn, a, qn, qn, an + 1, t + an + 1 - qn);
    else
        _flint_mpn_mulhigh_known_low(t, b, bn, q, qn, a, qn, qn, an + 1, t + an + 1 - qn);

    res = (t[an - qn] == 0) && flint_mpn_equal_p(t, a + qn, an - qn);

    TMP_END;
    return res;
}
