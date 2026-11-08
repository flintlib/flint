/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"
#if FLINT_BITS == 64
#include "mpn_extras/recip_small.h"
#endif

/*
    Arctangents of floating-point mantissas (mp_real.h), the kernels of
    atan and atan2 of nfloat, arf and arb.

    atan |x| for |x| < 1 is the kernel _mp_real_atan_kernel on |x| at
    enough limbs for the relative accuracy (the leading zero bits of |x|
    on top of nk); for |x| > 1, pi/2 - atan(1/|x|) with one approximate
    reciprocal at nk limbs.

    atan2(y, x) is c pi/4 + sigma atan(t) for the ratio t = min(|x|, |y|)
    / max(|x|, |y|) < 1 of the mantissas (one approximate division, with
    the truncation of the mantissas accounted for) and the octant c in
    {0, 2, 4}, sigma = +-1 (the caller sorts the arguments and applies
    the sign of y).  Only theta = atan(t) itself can be small (c = 0);
    then the kernel runs at as many more bits as t has leading zeros, and
    a tiny t gives theta = t (1 - delta) as the ratio itself.
*/

/* (y, w + 1) = c floor(pi/4 B^w) + sigma atan(v) for v at w fraction
   limbs within ev ulps (v < 1), c in {0, ..., 4}, sigma = +-1, no
   cancellation (c pi/4 >= 2 atan(v) when sigma < 0); atan' <= 1, and
   pi/4 is within one ulp below, c times */
static void
_finish(nn_ptr y, ulong * err, nn_srcptr v, slong w, ulong ev, int c, int sigma)
{
    ulong e;

    if (c == 0)
    {
        _mp_real_atan_kernel(y, &e, v, w);
        y[w] = 0;
        *err = e + ev;
    }
    else
    {
        nn_ptr t;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * w);
        _mp_real_atan_kernel(t, &e, v, w);

        y[w] = mpn_mul_1(y, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
        if (sigma > 0)
            y[w] += mpn_add_n(y, y, t, w);
        else
            y[w] -= mpn_sub_n(y, y, t, w);

        *err = e + ev + c;
        TMP_END;
    }
}

void
_mp_real_atan_fp(nn_ptr y, slong * wout, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax)
{
    slong nk = _mp_real_elem_limbs(pb);
    nn_ptr v;
    slong w;
    ulong ev;
    TMP_INIT;

    FLINT_ASSERT(wmax >= nk);

    if (e <= 0)
    {
        /* |x| < 1: the kernel on |x| at enough limbs for the leading
           zero bits (within the caller's buffers) */
        int trunc;

        w = FLINT_MIN(_mp_real_elem_limbs(pb - e), wmax);
        TMP_START;
        v = TMP_ALLOC(sizeof(ulong) * w);
        trunc = _mp_real_fp_get_fixed(v, w, w, d, n, e);
        _finish(y, err, v, w, trunc, 0, 1);
        TMP_END;
        *wout = w;
        return;
    }

    w = nk;
    *wout = w;

    if (e == 1 && d[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(d, n - 1))
    {
        /* |x| = 1: pi/4 */
        flint_mpn_copyi(y, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w);
        y[w] = 0;
        *err = 1;
        return;
    }

    TMP_START;
    v = TMP_ALLOC(sizeof(ulong) * (2 * w + 4));

    if (e > FLINT_BITS * (w + 1))
    {
        /* 1/|x| below one ulp: pi/2 - v with 0 < v < 2^(1-e) */
        _flint_mpn_zero_small(v, w);
        ev = 1;
    }
    else
    {
        slong ls, k;
        unsigned int bs;

        /* a mantissa beyond w + 1 limbs truncated: 1/D' - 1/D < 1/D'
           relatively, below 2 B^-(w + 1), times v <= 1: under one ulp */
        ev = 0;
        if (n > w + 1)
        {
            d += n - (w + 1);
            n = w + 1;
            ev = 1;
        }

#if FLINT_BITS == 64
        if (w >= 4 && w <= FLINT_MPN_RECIP_SMALL_MAX && w == n + 1)
        {
            /* 1/|x| = 2^(1 - e) / X, X = 2M = 1 + F in [1, 2) for the
               mantissa M: the dense reciprocal (I <= 1/X within six
               ulps), shifted right by e - 1 bits: v within seven ulps,
               below 1 */
            ulong F[FLINT_MPN_RECIP_SMALL_MAX], I[FLINT_MPN_RECIP_SMALL_MAX];
            slong s = e - 1;

            F[0] = 0;
            mpn_lshift(F + 1, d, n, 1);
            _flint_mpn_recip_small(I, F, w, FLINT_BITS * w - 3);
            ls = s / FLINT_BITS;
            bs = s % FLINT_BITS;
            for (k = 0; k < w; k++)
                v[k] = 0;
            if (ls < w)
            {
                if (bs != 0)
                    mpn_rshift(v, I + ls, w - ls, bs);
                else
                    for (k = 0; k < w - ls; k++)
                        v[k] = I[ls + k];
            }
            ev += 7;
            _finish(y, err, v, w, ev, 2, -1);
            TMP_END;
            return;
        }
#endif

        {
            /* v = 1/|x| = B^w 2^-e / M < 1: R = floor(B^(w + n) / D) or
               one more, in (B^w, 2 B^w], shifted right by e bits (the
               kernel costs the same for all v < 1, so this beats a
               reduction to (|x| - 1) / (|x| + 1) on [1, 2)) */
            nn_ptr R = v + w + 2;

            flint_mpn_invapprox(R, d, n, w + n);
            /* (R, w + 2), the top limb zero */
            ls = e / FLINT_BITS;
            bs = e % FLINT_BITS;
            ev += 2;
            if (ls >= w + 1)
                _flint_mpn_zero_small(v, w + 1);
            else
            {
                if (bs != 0)
                    mpn_rshift(R + ls, R + ls, w + 2 - ls, bs);
                flint_mpn_copyi(v, R + ls, w + 1 - ls);
                _flint_mpn_zero_small(v + w + 1 - ls, ls);
            }

            if (v[w] != 0)
            {
                /* v = 1 by the rounding of R (e = 1): one ulp below */
                FLINT_ASSERT(e == 1);
                for (k = 0; k < w; k++)
                    v[k] = UWORD_MAX;
                ev++;
            }
        }
    }

    _finish(y, err, v, w, ev, 2, -1);
    TMP_END;
}

/* (q, w + 1) = (A / D) B^w for the fractions A = (a, an), D = (d, dn) in
   [1/2, 1) (the top limbs set), the mantissas truncated to w limbs:
   returns the error bound in units.  |A/D - A'/D'| <= |A - A'| / D'
   + (A'/D') |D - D'| / D with D, D' >= 1/2 and A'/D' < 2, plus one for
   the floor and one for the approximate quotient. */
static ulong
_ratio(nn_ptr q, nn_srcptr a, slong an, nn_srcptr d, slong dn, slong w)
{
    ulong ea = 0, ed = 0;

    if (an > w)
    {
        a += an - w;
        an = w;
        ea = 1;
    }
    if (dn > w)
    {
        d += dn - w;
        dn = w;
        ed = 1;
    }

    flint_mpn_divapprox_fraction(q, a, an, d, dn, w - an + dn);
    return 2 * ea + 5 * ed + 2;
}

void
_mp_real_atan2_fp(nn_ptr y, slong * len, slong * ue, ulong * elo, ulong * ehi,
    nn_srcptr num, slong nn, nn_srcptr den, slong dn, slong d, int c, int sigma, slong pb, slong wmax)
{
    nn_ptr q, v;
    slong w, z, i, nk = _mp_real_elem_limbs(pb);
    ulong eq, ev, err;
    TMP_INIT;

    FLINT_ASSERT(d >= 0 && wmax >= nk + 1);

    /* the leading zero bits of t: d - 1 or d */
    z = (d >= 1) ? d - 1 : 0;

    if (c == 0)
    {
        /* tiny theta = t (1 - delta), delta < t^2 / 3 < 2^(-2z) */
        if (2 * z >= FLINT_BITS * (nk + 1) + 4)
        {
            /* the ratio of the mantissas at nk + 1 limbs, exactly floored
               (an approximate quotient would leave the side of a dyadic
               t unknown, which a correctly rounding caller needs):
               Q = floor((Mnum / Mden) B^w), t = Q 2^(-d - FLINT_BITS w)
               when the remainder vanishes, else t in (Q, Q + 1) units;
               delta t is below one unit */
            nn_ptr a, r;
            slong f, drop = 0;
            int exact;

            w = nk + 1;
            f = w - nn + dn;
            if (f < 0)
            {
                /* floor(floor(Mnum / B^-f) / Mden) */
                drop = -f;
                f = 0;
            }

            TMP_START;
            a = TMP_ALLOC(sizeof(ulong) * (nn - drop + f + dn));
            r = a + nn - drop + f;
            _flint_mpn_zero_small(a, f);
            flint_mpn_copyi(a + f, num + drop, nn - drop);
            flint_mpn_tdiv_qr(y, r, a, nn - drop + f, den, dn);
            exact = flint_mpn_zero_p(r, dn) && flint_mpn_zero_p(num, drop);
            TMP_END;

            *len = w + 1;
            *ue = -d - FLINT_BITS * w;
            *elo = 1;
            *ehi = exact ? 0 : 1;
            return;
        }

        w = FLINT_MIN(_mp_real_elem_limbs(pb + z), wmax);
    }
    else
    {
        w = nk;

        /* t below one ulp: theta = c pi/4 + sigma t with t < 2^-z, that
           is c pi/4 within one unit (as _finish with t = 0, whose kernel
           value and error vanish) */
        if (z >= FLINT_BITS * w + 2)
        {
            y[w] = mpn_mul_1(y, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
            err = 1 + c;
            goto done;
        }
    }

    TMP_START;
    q = TMP_ALLOC(sizeof(ulong) * (2 * w + 2));
    v = q + w + 1;

    /* the mantissa ratio at w limbs: Q = (Mnum / Mden) B^w in
       (B^w / 2, 2 B^w), t = Q B^-w 2^-d < 1 as a fraction at w limbs,
       shifted right by d bits (one more ulp for the shifted-out bits,
       the error scaled down) */
    eq = _ratio(q, num, nn, den, dn, w);

    {
        slong ls = d / FLINT_BITS;
        unsigned int bs = d % FLINT_BITS;

        _flint_mpn_zero_small(v, w + 1);
        for (i = 0; i <= w && ls + i <= w; i++)
        {
            v[i] = q[ls + i] >> bs;
            if (bs != 0 && ls + i + 1 <= w)
                v[i] |= q[ls + i + 1] << (FLINT_BITS - bs);
        }

        ev = (d >= FLINT_BITS) ? 2 : (eq >> d) + 2;
    }

    /* t < 1 exactly; the approximation may reach 1 */
    if (v[w] != 0)
    {
        for (i = 0; i < w; i++)
            v[i] = UWORD_MAX;
        ev++;
    }

    _finish(y, &err, v, w, ev, c, sigma);
    TMP_END;

done:
    *len = w + 1;
    *ue = -FLINT_BITS * w;
    *elo = *ehi = err;
}
