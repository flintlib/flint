/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* log of a positive ball to a relative accuracy of about 2^-prec.

   The midpoint m is evaluated exactly as given; the radius rho enters
   at the end as rho / (m - rho), which bounds |log y - log m| over the
   ball (mean value theorem), and prec is first lowered to what that
   allows.  With m = 2^E u, u in [1/2, 1):

   - near 1 (E = 0 or 1, D = m - 1 exact, |D| < 2^-z): D, D - D^2/2 or
     D - D^2/2 + D^3/3 with the remainder in the radius; from 32
     leading zero bits (by size) 2 atanh(D / (2 + D)) by the series of
     _mp_real_atanh_rs, after one ball division;
   - otherwise log m = (E - 1) log 2 + log1p(2u - 1) by the bitwise
     kernel up to 600 limbs, E log 2 - (-log u) by the Newton kernel
     beyond, the sum in ball arithmetic.  Near 1 either form cancels
     up to z + 1 bits (log 2 against log1p(2u - 1) for E = 0, against
     -log u for E = 1), which the working precision absorbs: the kernel
     runs at prec + z + 1 bits (absolute), so that the result keeps its
     relative accuracy. */

#define LG_BITWISE_MAX 600

/* the ratio and the atanh series (series_tapered.c / series_rs.c below
   32 zero bits, _mp_real_atanh_rs from 32) from z leading zero bits of
   D on, for m above 1 (E = 1) resp. below (E = 0), where the kernel
   (log1p of 2u - 1 close to 1) costs more (measured end to end, kernel
   / series time ratios: above 1, 1.2-1.4 at 5-6 limbs and 1.1-1.45 at
   10-24 limbs from z = 32, 1.1-1.2 at 7-8 limbs from z = 20, 1.15-1.25
   at 32-48 limbs from z = 64, 1.05-1.2 at 64-256 limbs from z = 96-256;
   below 1, 1.05-1.4 at 1-24 limbs from z = 14-28, 1.1-1.2 at 32-128
   limbs from z = 40-96, 1.05-1.25 at 192-512 limbs from z = 192-256) */
static slong
_lg_series_min_z(slong n, int E)
{
    if (E == 1)
    {
        if (n <= 4)
            return WORD_MAX;
        if (n <= 6)
            return 32;
        if (n <= 8)
            return 20;
        if (n <= 24)
            return 32;
        if (n <= 48)
            return 64;
        if (n <= 64)
            return 96;
        if (n <= 192)
            return 192;
        if (n <= 256)
            return 256;
        return WORD_MAX;
    }
    if (n <= 1)
        return 16;
    if (n <= 2)
        return 28;
    if (n <= 6)
        return 24;
    if (n <= 8)
        return 16;
    if (n <= 12)
        return 24;
    if (n <= 24)
        return 28;
    if (n <= 40)
        return 40;
    if (n <= 64)
        return 48;
    if (n <= 128)
        return 96;
    if (n <= 384)
        return 192;
    return 256;
}

/* kernel limbs for an absolute accuracy of about 2^-p: the bitwise
   bound 3 r + 64, the truncations and two bits */
static slong
_lg_limbs(slong p)
{
    slong n = (p + 12 + FLINT_BITS - 1) / FLINT_BITS, g;

    if (n <= LG_BITWISE_MAX)
        g = FLINT_BIT_COUNT(3 * (ulong) _mp_real_log1p_default_r_inline(n) + 72) + 3;
    else
        g = 8;
    n = (p + g + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* D = m - 1 exactly (ball arithmetic) */
static void
_lg_delta(mp_real_t D, const mp_real_t m)
{
    mp_real_t one;
    slong top = FLINT_MAX(m->exp, 1), bot = FLINT_MIN(m->exp - m->size, 0);

    mp_real_init(one);
    mp_real_set_ui(one, 1);
    mp_real_sub(D, m, one, top - bot + 1);
    FLINT_ASSERT(D->err == 0);
    mp_real_clear(one);
}

/* for m in [1/2, 2) (E = 0 or 1): z with |m - 1| < 2^-z, read off the
   mantissa (the leading zeros of m - 1's fraction for m >= 1, one less
   than the leading ones of m for m < 1); WORD_MAX for m = 1 */
static slong
_lg_closeness(const mp_real_t m, int E)
{
    slong j, z;

    if (E == 1)
    {
        /* top limb 1, the fraction below */
        for (j = m->size - 2; j >= 0 && m->d[j] == 0; j--)
            ;
        if (j < 0)
            return WORD_MAX;
        return FLINT_BITS * (m->size - 2 - j) + flint_clz(m->d[j]);
    }
    else
    {
        /* m = 0.11...10... : 1 - m <= 2^-(ones) */
        z = 0;
        for (j = m->size - 1; j >= 0 && m->d[j] == ~UWORD(0); j--)
            z += FLINT_BITS;
        if (j >= 0)
            z += flint_clz(~m->d[j]);
        return z - 1;
    }
}

/* s = |m - 1| as the fraction (s, n), truncated towards zero for m > 1
   and upwards (from the truncated m) for m < 1: within one ulp */
static void
_lg_distance(nn_ptr s, slong n, const mp_real_t m, int E)
{
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC((n + 1) * sizeof(ulong));
    _mp_real_elem_copy(t, n + 1, n, m);
    if (E == 1)
        flint_mpn_copyi(s, t, n);           /* drop the units limb */
    else
        mpn_neg(s, t, n);                   /* 1 - m, m < 1 */
    TMP_END;
}

/* u = m 2^-E in [1/2, 1) as the fraction (v, n), truncated, or with
   sh1 = 1 the fraction 2u - 1 (u's bits below the leading one, shifted
   up by one more); returns 1 if bits were dropped */
static int
_lg_mantissa(nn_ptr v, slong n, const mp_real_t m, int sh1)
{
    slong j, idx, size = m->size, base = size - n;
    int t = flint_clz(m->d[size - 1]) + sh1;
    ulong hi, lo;

    if (t == FLINT_BITS)
    {
        t = 0;
        base--;
    }

    /* v[j] takes mantissa limb base + j shifted left by t */
    for (j = 0; j < n; j++)
    {
        idx = base + j;
        hi = (idx >= 0) ? m->d[idx] : 0;
        lo = (idx >= 1) ? m->d[idx - 1] : 0;
        v[j] = (t == 0) ? hi : ((hi << t) | (lo >> (FLINT_BITS - t)));
    }

    /* the dropped bits: those of limb base - 1 not shifted in, and
       everything below */
    idx = base - 1;
    if (idx < 0)
        return 0;
    if (t != 0 && (m->d[idx] << t) != 0)
        return 1;
    if (t == 0 && m->d[idx] != 0)
        return 1;
    for (j = 0; j < idx; j++)
        if (m->d[j] != 0)
            return 1;
    return 0;
}

/* log(m) for exact m > 0, m != 1, to about prec bits */
static void
_lg_mid(mp_real_t res, const mp_real_t m, slong prec)
{
    slong E, z = 0, n;
    int near;

    E = FLINT_BITS * (m->exp - 1) + FLINT_BIT_COUNT(m->d[m->size - 1]);
    near = (E == 0 || E == 1);

    if (near)
    {
        z = _lg_closeness(m, (int) E);
        if (z == WORD_MAX)
        {
            /* m = 1 (the midpoint of an inexact ball may carry zero
               padding limbs) */
            mp_real_zero(res);
            return;
        }
    }

    if (near && z >= prec + 3)
    {
        /* |log(1 + D) - D| <= D^2 / (2 (1 - |D|)) < 2^(-2z) */
        mp_real_t D, zero;
        mp_real_init(D);
        mp_real_init(zero);
        _lg_delta(D, m);
        mp_real_add(res, D, zero, mp_real_prec_bits(prec + 4));
        mp_real_add_error_2exp_si(res, -2 * z);
        mp_real_clear(D);
        mp_real_clear(zero);
        return;
    }

    if (near && z >= 2 && 3 * z >= prec + 3)
    {
        /* log(1 + D) = D - D^2/2 (+ D^3/3), |D| = s < 2^-z, in fixed
           point at prec + z + 8 bits; the first omitted term and the
           rest, sum_{k > d} s^k / k < s^(d+1), in the radius (relative
           2^-(dz - 1) <= 2^-(prec + 2)); the truncations of s and of
           the powers stay below 5 ulps */
        int d = (2 * z >= prec + 3) ? 2 : 3;
        nn_ptr s, q, c;
        int neg = (E == 0);
        TMP_INIT;

        n = (prec + z + 8 + FLINT_BITS - 1) / FLINT_BITS;
        TMP_START;
        s = TMP_ALLOC(3 * n * sizeof(ulong));
        q = s + n;
        c = q + n;
        _lg_distance(s, n, m, (int) E);
        flint_mpn_sqrhigh(q, s, n);
        if (d == 3)
        {
            flint_mpn_mulhigh_n(c, q, s, n);
            mpn_divrem_1(c, 0, c, n, 3);             /* s^3/3 */
        }
        mpn_rshift(q, q, n, 1);                      /* s^2/2 */

        /* |log(1 + D)| = s - s^2/2 + s^3/3 (D > 0), s + s^2/2 + s^3/3
           (D < 0) */
        mp_real_fit_length(res, n + 1);
        flint_mpn_copyi(res->d, s, n);
        res->d[n] = 0;
        if (neg)
            mpn_add_n(res->d, res->d, q, n);
        else
            mpn_sub_n(res->d, res->d, q, n);
        if (d == 3)
            mpn_add_n(res->d, res->d, c, n);
        _mp_real_elem_finish(res, n, 5, neg);
        mp_real_add_error_2exp_si(res, -(d + 1) * z);
        TMP_END;
        return;
    }

    n = _lg_limbs(prec + (near ? z + 1 : 0));

    if (near && z >= _lg_series_min_z(n, (int) E))
    {
        /* log(1 + D) = 2 atanh(D / (2 + D)): w = s / (2 +- s) < 2^-z,
           by the series of _mp_real_atanh_rs (z >= 32) or the tapered
           series */
        nn_ptr s, w;
        ulong ea;
        int neg = (E == 0);
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(2 * n * sizeof(ulong));
        w = s + n;
        _lg_distance(s, n, m, (int) E);
        _mp_real_elem_half_ratio(w, s, n, neg);
        mp_real_fit_length(res, n + 1);
        if (z >= 32)
            _mp_real_atanh_rs(res->d, &ea, w, n);
        else
        {
            /* w < 2^-z */
            _mp_real_series_atan(res->d, w, n, (flint_bitcnt_t) z, MP_REAL_SERIES_ATANH);
            ea = 6;
        }
        res->d[n] = 0;
        /* s within 1 ulp, w within 2 (atanh' < 1 + 2^-63), doubled */
        _mp_real_elem_finish(res, n, ea + 3, neg);
        mp_real_mul_2exp_si(res, res, 1);
        TMP_END;
        return;
    }

    {
        /* log m = c log 2 + sigma F in fixed point: F = log1p(2u - 1)
           (sigma = 1, c = E - 1) or -log u (sigma = -1, c = E), and
           |c| L with L = floor(log 2 B^(n+1)) (read in place), whose error times
           |c| < 2^62 stays below a quarter ulp of B^-n */
        nn_ptr v, T;
        ulong err, cc, cy;
        slong c;
        int trunc, sigma, neg;
        TMP_INIT;

        TMP_START;
        v = TMP_ALLOC((2 * n + 2) * sizeof(ulong));
        T = v + n;

        /* 2u - 1 (bitwise) resp. u (Newton) */
        trunc = _lg_mantissa(v, n, m, n <= LG_BITWISE_MAX);

        /* c = 0 (E = 1 bitwise, E = 0 Newton): the kernel output is
           the result, written straight into res */
        if ((n <= LG_BITWISE_MAX) ? (E == 1) : (E == 0))
        {
            mp_real_fit_length(res, n + 1);
            if (n <= LG_BITWISE_MAX)
            {
                if (!_mp_real_log1p_opt(res->d, &err, v, n))
                    _mp_real_log1p_bitwise_rs(res->d, &err, v, n, 0);
            }
            else
                _mp_real_neglog_newton(res->d, &err, v, n);
            res->d[n] = 0;
            _mp_real_elem_finish(res, n, err + 2 * trunc, n > LG_BITWISE_MAX);
            TMP_END;
            return;
        }

        /* F, written one limb up into T = F B^-n at B^-(n+1) */
        if (n <= LG_BITWISE_MAX)
        {
            /* 2u - 1 is u's bits below the leading one (twice u's
               truncation error), log1p' <= 1 */
            if (!_mp_real_log1p_opt(T + 1, &err, v, n))
                _mp_real_log1p_bitwise_rs(T + 1, &err, v, n, 0);
            sigma = 1;
            c = E - 1;
        }
        else
        {
            /* (-log)' <= 2 on [1/2, 1) */
            _mp_real_neglog_newton(T + 1, &err, v, n);
            sigma = -1;
            c = E;
        }
        err += 2 * trunc;
        T[0] = 0;
        T[n + 1] = 0;

        /* T = |c| L +- F at B^-(n+1), n + 2 limbs, by one multiply-add
           resp. -subtract of L (n + 1 limbs, read in place) */
        cc = (c < 0) ? -(ulong) c : (ulong) c;
        {
            nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, n + 1);

            if ((sigma > 0 && c >= 0) || (sigma < 0 && c <= 0))
            {
                /* |c| L + F; negative for sigma = -1 (c <= 0) */
                MP_REAL_ADDMUL_1(cy, T, L, n + 1, cc);
                T[n + 1] += cy;
                neg = (sigma < 0);
            }
            else
            {
                /* F - |c| L < 0 in two's complement, negated: |c| L - F;
                   negative for sigma = 1 (c < 0) */
                MP_REAL_SUBMUL_1(cy, T, L, n + 1, cc);
                T[n + 1] -= cy;
                mpn_neg(T, T, n + 2);
                neg = (sigma > 0);
            }
        }

        /* drop the extra fraction limb (one more ulp) */
        mp_real_fit_length(res, n + 1);
        flint_mpn_copyi(res->d, T + 1, n + 1);
        _mp_real_elem_finish(res, n, err + 2, neg);
        TMP_END;
    }
}

int
mp_real_log_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    ulong xerr = x->err;
    slong xanc = x->exp - x->size;

    prec = FLINT_MAX(prec, 2);

    /* strictly positive: the midpoint's mantissa exceeds the radius
       count at the common anchor */
    if (x->size == 0 || x->negative
        || (xerr != 0 && x->size == 1 && x->d[0] <= xerr))
    {
        mp_real_zero(res);
        return 0;
    }

    /* log 1 = 0 */
    if (xerr == 0 && x->size == 1 && x->exp == 1 && x->d[0] == 1)
    {
        mp_real_zero(res);
        return 1;
    }

    mid = *x;
    mid.err = 0;

    if (xerr == 0)
    {
        _lg_mid(res, &mid, prec);
        return 1;
    }
    else
    {
        /* |log y - log m| <= rho / (m - rho) < 2^(rel + 1) for
           rho < 2^rel m, rel <= -1; relative to |log m| >= 2^lg */
        slong rel = mp_real_rel_radius_lt_2exp_si(x), lg, acc, E;

        E = FLINT_BITS * (x->exp - 1) + FLINT_BIT_COUNT(x->d[x->size - 1]);
        if (E == 0 || E == 1)
        {
            slong z = _lg_closeness(&mid, (int) E);
            lg = (z == WORD_MAX) ? -FLINT_BITS * x->size : -z - 2;
        }
        else
            lg = FLINT_BIT_COUNT(FLINT_ABS(E) - 1) - 2;

        acc = -(FLINT_MAX(rel, -WORD_MAX / 4) + 1) - lg;
        prec = FLINT_MIN(prec, FLINT_MAX(acc + 8, 10));

        /* rho / (m - rho) = (rd / (mlo (1 - r))) B^(xanc - exp + 1) with
           m >= mlo B^(exp - 1) and r >= rho / (mlo B^(exp - 1)), in
           double arithmetic rounded up; read before res is written (it
           may alias x) */
        {
            double rd = (double) xerr, mlo = _mp_real_mag_lo(x);
            slong a = xanc - (x->exp - 1);
            /* rd / mlo in [2^-64, 2^64]: from a = 2 on, r >= 2^64 (the
               ball is wide); below, scaled without over- or underflow */
            double r = (a >= 2) ? 0x1p64
                : _mp_real_elem_scale_up(rd / mlo, a) * (1.0 + 0x1p-50);

            if (r <= 0.5)
            {
                double v = rd / (mlo * (1.0 - r)) * (1.0 + 0x1p-48);
                _lg_mid(res, &mid, prec);
                _mp_real_elem_add_rad_d(res, v, a);
            }
            else
            {
                /* a wide ball: the bound in ball arithmetic */
                mp_real_t R, lo, Q;
                mp_real_init(R);
                mp_real_init(lo);
                mp_real_init(Q);
                _mp_real_set_mpn_2exp(R, &xerr, 1, FLINT_BITS * xanc);
                mp_real_sub(lo, &mid, R, x->size + 2);
                FLINT_ASSERT(lo->err == 0 && !lo->negative && lo->size != 0);
                mp_real_div(Q, R, lo, 2);
                _lg_mid(res, &mid, prec);
                _mp_real_elem_add_mag(res, Q);
                mp_real_clear(R);
                mp_real_clear(lo);
                mp_real_clear(Q);
            }
        }
        return 1;
    }
}
