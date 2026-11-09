/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"

/*
    Logarithms of floating-point mantissas (mp_real.h), the kernels of
    the logarithms of nfloat, arf and arb.

    For u = M 2^E, M in [1/2, 1), log u = (E - 1) log 2 + log1p(2M - 1),
    with the log1p kernel on v = 2M - 1 in [0, 1) at w fraction limbs
    and log 2 the floor L of log 2 B^(w + 1), whose error times
    |E - 1| < 2^(FLINT_BITS - 3) stays below an eighth of an ulp of
    B^-w (a bignum E - 1 takes L at as many more limbs).  Near 1 (E = 0
    or 1) with |u - 1| < 2^-z the sum cancels up to z + 1 bits, which the
    kernel absorbs by running at z + 1 more bits; from enough leading
    zero bits (by kernel limbs, the crossover of log.c),
    log(1 +- s) = +-2 atanh(s / (2 +- s)) by the atanh series instead.

    log2 u = (E - 1) + log1p(2M - 1) / log 2, the division a
    multiplication by 1/log 2 (the tabulated constant), exact in its
    integral part.

    log1p(x) for |x| < 1/2 goes to the log(1 +- s) kernel on s = |x|;
    for other x the sum 1 + x is formed exactly as a mantissa.
*/

/* T = log2(e) * (Y, w) = Y + Y R at the same scale, R the fraction of
   1/log 2: (T, w + 1) with err ulps updated (the value 1.443.. times
   the error, plus the truncations) */
static void
_mul_ilog2(nn_ptr T, nn_srcptr Y, slong w, ulong * err)
{
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC((w + 1) * sizeof(ulong));
    flint_mpn_mulhigh_n(t, Y, _mp_real_const_ptr(MP_REAL_CONST_ID_INV_LOG2_FRAC, w), w);
    T[w] = mpn_add_n(T, Y, t, w);
    *err = *err + (*err >> 1) + 3;
    TMP_END;
}

/*
    log(1 + s) (neg = 0) resp. log(1 - s) (neg = 1) for the fraction
    (s, w), s < 1 resp. s <= 1/2, s known within serr ulps (as a lower
    bound for neg = 0, an upper bound for neg = 1), w at least the kernel
    limbs for z + 1 extra bits where s < 2^-z. Sets (y, w + 1), returns the
    exponent of its unit and sets *err.
*/
static slong
_log1pm_fixed(nn_ptr y, ulong * err, nn_srcptr s, ulong serr, slong w, int neg)
{
    slong z = _mp_real_elem_lzb(s, w);
    TMP_INIT;

    FLINT_ASSERT(z != WORD_MAX);

    if (z >= _mp_real_log_series_min_z(w, !neg) && z >= 8)
    {
        /* +-2 atanh(t), t = s / (2 +- s) within serr + 1 ulps, atanh' < 1 +
           2^-15 */
        nn_ptr t;
        ulong ea;

        TMP_START;
        t = TMP_ALLOC(w * sizeof(ulong));
        _mp_real_elem_half_ratio(t, s, w, neg);
        if (z >= 32)
            _mp_real_atanh_rs(y, &ea, t, w);
        else
        {
            _mp_real_series_atan(y, t, w, (flint_bitcnt_t) z, MP_REAL_SERIES_ATANH);
            ea = 6;
        }
        y[w] = 0;
        *err = ea + serr + 2;
        TMP_END;
        return 1 - FLINT_BITS * w;
    }
    else if (!neg)
    {
        /* the kernel (log1p' <= 1) */
        _mp_real_log1p_kernel(y, err, s, w);
        y[w] = 0;
        *err += serr;
        return -FLINT_BITS * w;
    }
    else
    {
        /* -log(1 - s) = log 2 - log1p(1 - 2s), 1 - 2s within 2 serr ulps,
           log1p' <= 1; F at w fraction limbs, L at w + 1 */
        nn_ptr v, T;
        ulong ek;

        TMP_START;
        v = TMP_ALLOC((2 * w + 3) * sizeof(ulong));
        T = v + w + 1;

        v[w] = mpn_lshift(v, s, w, 1);
        mpn_neg(v, v, w);

        T[0] = 0;
        _mp_real_log1p_kernel(T + 1, &ek, v, w);
        mpn_sub_n(T, _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, w + 1), T, w + 1);
        _flint_mpn_copy_small(y, T + 1, w);
        y[w] = 0;
        /* L within one ulp of B^(w+1), the dropped limb one ulp */
        *err = ek + 2 * serr + 2;
        TMP_END;
        return -FLINT_BITS * w;
    }
}

/* the kernel limbs with extra more bits (at least two limbs on 32-bit
   machines, as the kernels require) */
static slong
_limbs(slong pb, slong extra)
{
    return _mp_real_elem_limbs(pb + extra);
}

int
_mp_real_log_fp(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, slong e, int base2, slong pb)
{
    nn_ptr s;
    slong w, z;
    int trunc, neg;
    TMP_INIT;

    /* (|e| <= 2^(FLINT_BITS - 3), NFLOAT_MIN_EXP included: |c| <=
       2^(FLINT_BITS - 3) + 1 keeps |c| L within the integral limb and
       |c| (log 2 - L) below a unit) */
    FLINT_ASSERT(FLINT_ABS(e) <= (WORD(1) << (FLINT_BITS - 3)));

    if (e == 0 || e == 1)
    {
        /* s = |u - 1| with u - 1 = (2M - 1) for E = 1, -(1 - M) for E = 0 */
        neg = (e == 0);

        if (!neg)
        {
            /* the fraction bits of 2M, exactly */
            nn_ptr t;
            slong j;

            TMP_START;
            t = TMP_ALLOC((n + 1) * sizeof(ulong));
            _mp_real_fp_get_fixed(t, n + 1, n, d, n, 1);
            for (j = n - 1; j >= 0 && t[j] == 0; j--)
                ;
            if (j < 0)
            {
                /* u = 1: log u = 0 exactly */
                *len = 0;
                *ue = 0;
                *err = 0;
                TMP_END;
                return 0;
            }
            z = FLINT_BITS * (n - 1 - j) + flint_clz(t[j]);
            TMP_END;
        }
        else
        {
            /* 1 - M: the leading ones of M */
            slong j;
            z = 0;
            for (j = n - 1; j >= 0 && d[j] == UWORD_MAX; j--)
                z += FLINT_BITS;
            if (j >= 0)
                z += flint_clz(~d[j]);
            z--;
        }

        /* s at w fraction limbs for z + 1 extra bits (the cancellation
           cannot exceed the bits of the mantissa) */
        w = _limbs(pb, FLINT_MIN(z, FLINT_BITS * n) + 1 + 2 * base2);

        TMP_START;
        s = TMP_ALLOC((w + 2) * sizeof(ulong));

        if (!neg)
        {
            trunc = _mp_real_fp_get_fixed(s, w + 1, w, d, n, 1);
            s[w] = 0;
        }
        else
        {
            /* 1 - M rounded up: the negated truncation of M (s[0] is
               written first only to keep gcc from warning) */
            s[0] = 0;
            trunc = _mp_real_fp_get_fixed(s, w, w, d, n, 0);
            mpn_neg(s, s, w);
        }

        *ue = _log1pm_fixed(y, err, s, trunc, w, neg);
        *len = w + 1;

        if (base2)
        {
            nn_ptr T = TMP_ALLOC((w + 3) * sizeof(ulong));
            _mul_ilog2(T, y, w + 1, err);
            _flint_mpn_copy_small(y, T, w + 2);
            *len = w + 2;
        }

        TMP_END;
        return neg;
    }
    else
    {
        /* (E - 1) log 2 + F, F = log1p(2M - 1), in fixed point with one
           integral limb and w + 1 fraction limbs: T = |c| L +- F */
        nn_ptr T;
        ulong cc, cy;
        slong c = e - 1;

        w = _limbs(pb, 0);

        TMP_START;
        s = TMP_ALLOC((2 * w + 6) * sizeof(ulong));
        T = s + w + 2;

        trunc = _mp_real_fp_get_fixed(s, w + 1, w, d, n, 1);   /* 2M at w limbs */
        T[0] = 0;
        _mp_real_log1p_kernel(T + 1, err, s, w);
        T[w + 1] = 0;
        *err += 2 * trunc;
        cc = (c < 0) ? -(ulong) c : (ulong) c;

        if (base2)
        {
            /* (E - 1) + F / log 2: the integer part is exact */
            _mul_ilog2(y, T + 1, w, err);
            /* G = F / log 2 < 1 (up to the error) at w fraction limbs */
            if (c > 0)
            {
                y[w] += cc;
                neg = 0;
            }
            else
            {
                /* |c| - G with |c| >= 2 */
                mpn_neg(y, y, w + 1);
                y[w] += cc;
                neg = 1;
            }
            *len = w + 1;
            *ue = -FLINT_BITS * w;
        }
        else
        {
            nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, w + 1);

            if (c > 0)
            {
                MP_REAL_ADDMUL_1(cy, T, L, w + 1, cc);
                T[w + 1] += cy;
                neg = 0;
            }
            else
            {
                /* F - |c| L < 0, negated */
                MP_REAL_SUBMUL_1(cy, T, L, w + 1, cc);
                T[w + 1] -= cy;
                mpn_neg(T, T, w + 2);
                neg = 1;
            }

            /* |c| (log 2 - L) < 2^(FLINT_BITS - 3) B^-(w+1) and the dropped limb */
            *err += 2;
            _flint_mpn_copy_small(y, T + 1, w + 1);
            *len = w + 1;
            *ue = -FLINT_BITS * w;
        }

        TMP_END;
        return neg;
    }
}

int
_mp_real_log_fp_huge(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, nn_srcptr c, slong cn, int cneg, slong pb)
{
    nn_ptr s, T, P;
    nn_srcptr L;
    slong w, wl;
    int trunc;
    TMP_INIT;

    /* |c| >= 2, so that |c| L - F >= 2 L - log 2 > 0 */
    FLINT_ASSERT(cn >= 1 && c[cn - 1] != 0);
    FLINT_ASSERT(cn > 1 || c[0] >= 2);

    /* F = log1p(2M - 1) at w fraction limbs; |c| L with L at wl = w + cn + 1
       limbs, so that |c| (log 2 - L) < B^cn B^-wl = B^-(w+1) */
    w = _limbs(pb, 0);
    wl = w + cn + 1;

    TMP_START;
    s = TMP_ALLOC((w + 2 + w + 2 + cn + wl) * sizeof(ulong));
    T = s + w + 2;
    P = T + w + 2;

    trunc = _mp_real_fp_get_fixed(s, w + 1, w, d, n, 1);
    _mp_real_log1p_kernel(T, err, s, w);
    *err += 2 * trunc;

    /* P = |c| L with cn integral and wl fraction limbs, +- F at the top
       w of the fraction limbs (no sign change: F < log 2 and |c| >= 2) */
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, wl);
    flint_mpn_mul(P, L, wl, c, cn);
    if (!cneg)
        mpn_add(P + wl - w, P + wl - w, cn + w, T, w);
    else
        mpn_sub(P + wl - w, P + wl - w, cn + w, T, w);

    /* the result at the unit B^-w: the dropped limbs below and
       |c| (log 2 - L) each below a unit */
    _flint_mpn_copy_small(y, P + wl - w, cn + w);
    *len = cn + w;
    *ue = -FLINT_BITS * w;
    *err += 2;
    TMP_END;
    return cneg;
}

int
_mp_real_log1p_fp(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb)
{
    TMP_INIT;

    FLINT_ASSERT(!(sgnbit && e >= 1));

    if (e <= -1)
    {
        /* |x| < 1/2: s = |x| directly */
        nn_ptr s;
        slong w;
        int trunc;

        w = _limbs(pb, -e + 1);
        TMP_START;
        s = TMP_ALLOC((w + 1) * sizeof(ulong));
        trunc = _mp_real_fp_get_fixed(s, w, w, d, n, e);
        *ue = _log1pm_fixed(y, err, s, trunc, w, sgnbit);
        *len = w + 1;
        TMP_END;
        return sgnbit;
    }

    if (sgnbit)
    {
        /* x in (-1, -1/2]: u = 1 - |x| in (0, 1/2] exactly */
        nn_ptr t;
        slong j, un;
        unsigned int c;
        int r;

        TMP_START;
        t = TMP_ALLOC(n * sizeof(ulong));
        mpn_neg(t, d, n);
        un = n;
        while (t[un - 1] == 0)
            un--;
        c = flint_clz(t[un - 1]);
        if (c != 0)
            mpn_lshift(t, t, un, c);
        /* exact values carry no low zero limbs */
        for (j = 0; t[j] == 0; j++)
            ;
        r = _mp_real_log_fp(y, len, ue, err, t + j, un - j, -(slong) (FLINT_BITS * (n - un) + c), 0, pb);
        TMP_END;
        return r;
    }
    else
    {
        slong w = _limbs(pb, 0);

        /* log1p(x) = log(x) + log1p(1/x) with 0 < log1p(1/x) < 2^(1-e),
           below one unit of the fixed-point logarithm */
        if (e > FLINT_BITS * w + 2)
        {
            int r = _mp_real_log_fp(y, len, ue, err, d, n, e, 0, pb);
            *err += 1;
            return r;
        }

        /* x >= 1/2: u = 1 + x exactly, as an integer with fn fraction
           limbs (x's unit is 2^(e - FLINT_BITS n), e >= 0) */
        {
            nn_ptr t;
            slong s = e - FLINT_BITS * n;
            slong fn = (s >= 0) ? 0 : (-s + FLINT_BITS - 1) / FLINT_BITS;
            slong tl = fn + e / FLINT_BITS + 2, j;
            unsigned int c;
            int r;

            TMP_START;
            t = TMP_ALLOC(tl * sizeof(ulong));
            _mp_real_fp_get_fixed(t, tl, fn, d, n, e);
            mpn_add_1(t + fn, t + fn, tl - fn, 1);
            while (t[tl - 1] == 0)
                tl--;
            c = flint_clz(t[tl - 1]);
            if (c != 0)
                mpn_lshift(t, t, tl, c);
            for (j = 0; t[j] == 0; j++)
                ;
            r = _mp_real_log_fp(y, len, ue, err, t + j, tl - j, FLINT_BITS * (tl - fn) - c, 0, pb);
            TMP_END;
            return r;
        }
    }
}
