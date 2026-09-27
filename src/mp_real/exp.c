/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"


/* exp of a ball to a relative accuracy of about 2^-prec.

   The midpoint m is evaluated exactly as given; the radius rho of x
   enters at the end as the factor [1 +- (e^rho - 1)], and prec is first
   lowered to what rho allows (the output's relative accuracy is about
   rho itself).  Small |m| take 1 + m (+ m^2/2), the series of
   _mp_real_exp_reduced (m > 0) or cosh - sinh (m < 0); otherwise

       |m| = q log 2 + t,   t in [0, log 2),

   with q a lower approximation of |m|/log 2 (two limbs of 1/log 2 and
   three umul_ppmm, at most one low, then corrected) and log 2 the
   cached floor at n + 1 fraction limbs; exp(m) = 2^q exp(t) for m > 0
   and 2^-(q+1) exp(log 2 - t) for m < 0, with the kernel on [0, 1).
   The reduction errs by q (log 2 - L) < 2^60 B^-(n+1) plus the
   truncations: under two ulps of B^-n in t, four in exp(t) < 2.

   The result exponent is about m / log 2 bits, so |m| is limited to
   2^(FLINT_BITS - 5) (a result exponent within the safe range);
   larger arguments throw. */

#define EX_BITWISE_MAX 600
#define EX_DIOPHANTINE_MAX 65536

/* the series of _mp_real_exp_reduced from z leading zero bits on, for
   m > 0 (measured end to end, kernel / series time ratios: 1.3-2 at
   2-10 limbs from z = 32, 1.1-1.6 at 12-40 limbs from z = 40-64,
   1.05-1.2 at 48-128 limbs from z = 96-192; below 32 zero bits the
   series never wins against the bitwise kernel) */
static slong
_ex_series_min_z(slong n)
{
    if (n <= 10)
        return 32;
    if (n <= 12)
        return 40;
    if (n <= 40)
        return 64;
    if (n <= 48)
        return 96;
    if (n <= 64)
        return 128;
    if (n <= 128)
        return 192;
    return WORD_MAX;
}

/* for m < 0, the alternating series of exp(-|m|) (series_rs.c; the
   hardcoded cosh - sinh at z >= 32 up to 8 limbs) from z on (measured
   end to end: 1.8-2.5 at 1-6 limbs from z = 32, 1.1-1.5 at 7-24 limbs
   from z = 22-28, 1.1-1.5 at 32-128 limbs from z = 40-96, 1.15-1.35 at
   192-256 limbs from z = 192) */
static slong
_ex_neg_series_min_z(slong n)
{
    if (n <= 6)
        return 32;
    if (n <= 8)
        return 22;
    if (n <= 24)
        return 28;
    if (n <= 40)
        return 40;
    if (n <= 48)
        return 48;
    if (n <= 128)
        return 96;
    if (n <= 256)
        return 192;
    return WORD_MAX;
}

/* (y, n + 1) = exp(v), v in [0, 1) at n fraction limbs */
static void
_ex_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n)
{
    slong z = _mp_real_elem_lzb(v, n);

    if (z == WORD_MAX)
    {
        flint_mpn_zero(y, n);
        y[n] = 1;
        *err = 0;
    }
    else if (z >= _ex_series_min_z(n))
        _mp_real_exp_reduced(y, err, v, n, (flint_bitcnt_t) z, 0);
    else if (n <= EX_BITWISE_MAX)
    {
        if (!_mp_real_exp_opt(y, err, v, n))
            _mp_real_exp_bitwise_rs(y, err, v, n, 0);
    }
    else if (n <= EX_DIOPHANTINE_MAX)
        _mp_real_exp_diophantine(y, err, v, n);
    else
        _mp_real_exp_notab(y, err, v, n);
}

/* kernel limbs for a relative accuracy of about 2^-p of exp(t) in
   [1, 2): the bitwise bound 9 r + 100 plus the reduction's 4 ulps and
   two bits */
static slong
_ex_limbs(slong p)
{
    slong n = (p + 14 + FLINT_BITS - 1) / FLINT_BITS, g;

    if (n <= EX_BITWISE_MAX)
        g = FLINT_BIT_COUNT(9 * (ulong) _mp_real_exp_default_r_inline(n) + 104) + 3;
    else
        g = 10;
    n = (p + g + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* R = floor((1/log 2 - 1) B^2), two limbs: 1/log 2 = 1.4426950408889634... */
#if FLINT_BITS == 64
#define EX_ILOG2_1 UWORD(0x71547652b82fe177)
#define EX_ILOG2_0 UWORD(0x7d0ffda0d23a7d11)
#else
#define EX_ILOG2_1 UWORD(0x71547652)
#define EX_ILOG2_0 UWORD(0xb82fe177)
#endif

/* X -= q L resp. X += q L over N limbs, returning the borrow resp.
   carry limb; unrolled in registers for constant N */
/* the limb of |m| at position i of the frame with the mantissa's limb 0
   at position sh */
FLINT_FORCE_INLINE ulong
_ex_frame_limb(const mp_real_t m, slong sh, slong i)
{
    i -= sh;
    return (i >= 0 && i < m->size) ? m->d[i] : 0;
}

/* exp(m) for exact m != 0 with |m| < 2^(FLINT_BITS - 5), into res at
   n kernel limbs.  X is |m| (m > 0) resp. its two's complement (m < 0)
   in a frame of N = n + 1 fraction limbs and one integral limb, the
   negation done in the copy; with q <= |m| / log 2 (exact or one low)
   one pass X - q L resp. X + (q + 1) L gives t = |m| - q L resp.
   t = (q + 1) L - |m| up to one more L, so that m = k L + t with
   t in [0, L], k = q resp. -(q + 1).  L = floor(log 2 B^N) is read in
   place (static or cached). */
static void
_ex_reduce(mp_real_t res, const mp_real_t m, slong n)
{
    slong N = n + 1, sh, q, k, j;
    nn_ptr X;
    nn_srcptr L;
    ulong x1, x0, h, l, c, u, cy, err;
    int neg = m->negative;
    TMP_INIT;

    FLINT_ASSERT(m->exp <= 1);

    TMP_START;
    X = TMP_ALLOC((N + 1) * sizeof(ulong));
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);
    sh = m->exp - m->size + N;

    /* q <= |m| / log 2: with X = x1 + x0/B and 1/log 2 = 1 + R/B^2,
       R = (R1, R0), the product is x1 + h1 + (x0 + l1 + hi(x1 R0))/B
       plus positive terms below 2/B (h1 B + l1 = x1 R1); dropping them
       and flooring gives the exact quotient or one less */
    x1 = _ex_frame_limb(m, sh, N);
    x0 = _ex_frame_limb(m, sh, N - 1);
    umul_ppmm(h, l, x1, EX_ILOG2_1);
    umul_ppmm(c, u, x1, EX_ILOG2_0);
    (void) u;
    {
        ulong s1, s0;
        add_ssaaaa(s1, s0, UWORD(0), x0, UWORD(0), l);
        add_ssaaaa(s1, s0, s1, s0, UWORD(0), c);
        (void) s0;
        q = (slong) (x1 + h + s1);
    }

    if (!neg)
    {
        /* X = |m| (truncated), t = X - q L */
        _mp_real_elem_copy(X, N + 1, N, m);
        MP_REAL_SUBMUL_1(cy, X, L, N, (ulong) q);
        X[N] -= cy;
        if (X[N] != 0 || mpn_cmp(X, L, N) >= 0)
        {
            /* q was one low */
            X[N] -= mpn_sub_n(X, X, L, N);
            q++;
        }
        k = q;
    }
    else
    {
        /* X = B^(N+1) - |m| (|m| truncated), t = X + (q + 1) L */
        slong lo = FLINT_MAX(sh, 0), drop = FLINT_MAX(-sh, 0);
        slong len = m->size - drop;

        flint_mpn_zero(X, lo);
        mpn_neg(X + lo, m->d + drop, len);      /* nonzero: borrows */
        for (j = lo + len; j <= N; j++)
            X[j] = ~UWORD(0);
        q++;
        MP_REAL_ADDMUL_1(cy, X, L, N, (ulong) q);
        X[N] += cy;
        if (X[N] != 0)
        {
            /* negative (q was one low): one more L */
            X[N] += mpn_add_n(X, X, L, N);
            q++;
        }
        k = -q;
    }
    FLINT_ASSERT(X[N] == 0);

    /* the kernel on the top n fraction limbs, straight into res; the
       reduction errs by (q + 1)(log 2 - L) < 2^60 B^-N, the truncations
       of m and of t: under two ulps of B^-n in t, four in exp(t) < 2 */
    mp_real_fit_length(res, n + 1);
    _ex_kernel(res->d, &err, X + 1, n);
    _mp_real_elem_finish(res, n, err + 4, 0);
    mp_real_mul_2exp_si(res, res, k);

    TMP_END;
}

/* exp(m), m exact, nonzero, |m| < 2^(FLINT_BITS - 5), to about prec
   bits */
static void
_ex_mid(mp_real_t res, const mp_real_t m, slong prec)
{
    slong emid, z, n;

    emid = FLINT_BITS * (m->exp - 1) + FLINT_BIT_COUNT(m->d[m->size - 1]);
    z = -emid;          /* |m| < 2^-z */

    if (z >= prec + FLINT_BITS)
    {
        /* |exp(m) - 1| <= |m| + m^2 < 2^(emid + 1) */
        _mp_real_elem_set_error(res, 1, emid + 1);
        return;
    }

    if (4 * z >= prec + 4)
    {
        /* 1 + m (+ m^2/2 (+ m^3/6)) in fixed point at n limbs, the
           first omitted term |m|^(d+1)/(d+1)! e^|m| < 2^((d+1) emid)
           in the radius; the truncations of m (twice, through exp' < 2)
           and of the powers stay below 6 ulps */
        int d = (2 * z >= prec + 4) ? 1 : (3 * z >= prec + 4) ? 2 : 3;
        nn_ptr v, a, q;
        int trunc;
        TMP_INIT;

        n = (prec + 12 + FLINT_BITS - 1) / FLINT_BITS;
        TMP_START;
        v = TMP_ALLOC(3 * n * sizeof(ulong));
        a = v + n;
        q = a + n;

        trunc = _mp_real_elem_copy(v, n, n, m);
        flint_mpn_copyi(a, v, n);                   /* A = v + v^3/6 */
        if (d >= 2)
        {
            flint_mpn_sqrhigh(q, v, n);
            mpn_rshift(q, q, n, 1);                 /* v^2/2 */
            if (d >= 3)
            {
                nn_ptr c = TMP_ALLOC(n * sizeof(ulong));
                flint_mpn_mulhigh_n(c, q, v, n);
                mpn_divrem_1(c, 0, c, n, 3);        /* v^3/6 */
                mpn_add_n(a, a, c, n);
            }
        }

        mp_real_fit_length(res, n + 1);
        flint_mpn_zero(res->d, n);
        res->d[n] = 1;
        if (m->negative)
            mpn_sub(res->d, res->d, n + 1, a, n);
        else
            mpn_add(res->d, res->d, n + 1, a, n);
        if (d >= 2)
            mpn_add(res->d, res->d, n + 1, q, n);
        _mp_real_elem_finish(res, n, 4 + 2 * trunc, 0);
        mp_real_add_error_2exp_si(res, (d + 1) * emid);
        TMP_END;
        return;
    }

    n = _ex_limbs(prec);

    if (m->exp <= 0 && !m->negative)
    {
        /* m in (0, 1): the kernel directly, m's limbs read in place
           when they can be */
        nn_srcptr v;
        nn_ptr buf;
        ulong err;
        int trunc;
        TMP_INIT;

        TMP_START;
        buf = TMP_ALLOC(n * sizeof(ulong));
        v = _mp_real_elem_frame(buf, n, m, res->d, NULL, &trunc);
        mp_real_fit_length(res, n + 1);
        _ex_kernel(res->d, &err, v, n);
        /* truncating m costs one ulp, exp(m) < e times that */
        _mp_real_elem_finish(res, n, err + 3 * trunc, 0);
        TMP_END;
        return;
    }

    if (m->exp <= 0 && z >= _ex_neg_series_min_z(n))
    {
        nn_ptr v, sh, ch;
        ulong err;
        int trunc;
        TMP_INIT;

        TMP_START;
        v = TMP_ALLOC((3 * n + 2) * sizeof(ulong));
        sh = v + n;
        ch = sh + n + 1;
        trunc = _mp_real_elem_copy(v, n, n, m);
        mp_real_fit_length(res, n + 1);
        if (z >= 32 && n <= 8)
        {
            /* m in (-2^-32, 0), small n: exp(m) = cosh |m| - sinh |m|
               by the hardcoded series */
            _mp_real_sinh_cosh_rs(sh, ch, &err, v, n);
            mpn_sub_n(res->d, ch, sh, n + 1);
            err *= 2;
        }
        else
        {
            /* the alternating series of exp(-|m|) */
            _mp_real_series_rs(res->d, v, n, (flint_bitcnt_t) z, MP_REAL_SERIES_EXP_NEG);
            err = 5;
        }
        _mp_real_elem_finish(res, n, err + 2 * trunc, 0);
        TMP_END;
        return;
    }

    _ex_reduce(res, m, n);
}

void
mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    slong e;
    ulong xerr = x->err;
    slong xanc = x->exp - x->size;

    prec = FLINT_MAX(prec, 2);

    if (x->size == 0 && x->err == 0)
    {
        mp_real_set_ui(res, 1);
        return;
    }

    /* |x| < 2^e over the whole ball */
    e = mp_real_abs_bound_lt_2exp_si(x);

    if (e > FLINT_BITS - 5)
        flint_throw(FLINT_ERROR, "mp_real_exp_bits: |x| >= 2^%wd, the result "
            "exponent would leave the safe range\n", (slong) (FLINT_BITS - 5));

    /* a ball around zero: exp in [1 +- 2^(e+1)] for |x| < 2^e <= 1/4
       (e^d - 1 <= 1.2 d, 1 - e^-d <= d), else [0 +- e^(2^e)] */
    if (x->size == 0)
    {
        if (e <= -2)
            _mp_real_elem_set_error(res, 1, e + 1);
        else
        {
            mp_real_t t, u;
            mp_real_init(t);
            mp_real_init(u);
            mp_real_set_ui(t, 1);
            mp_real_mul_2exp_si(t, t, e);
            mp_real_exp_bits(u, t, 30);
            mp_real_zero(res);
            _mp_real_elem_add_mag(res, u);
            mp_real_clear(t);
            mp_real_clear(u);
        }
        return;
    }

    mid = *x;
    mid.err = 0;

    if (xerr != 0)
    {
        /* rho < 2^(rel + emid); the output's relative accuracy is
           about rho */
        slong rel = mp_real_rel_radius_lt_2exp_si(x);
        slong emid = FLINT_BITS * (x->exp - 1) + FLINT_BIT_COUNT(x->d[x->size - 1]);
        slong acc = -(rel + emid);

        if (acc < 2)
        {
            /* exp over the ball lies in (0, exp(m + rho)] */
            mp_real_t hi, r, u;
            mp_real_init(hi);
            mp_real_init(r);
            mp_real_init(u);
            _mp_real_set_mpn_2exp(r, &xerr, 1, FLINT_BITS * xanc);
            mp_real_add(hi, &mid, r, x->size + 2);
            hi->err = 0;        /* an upper bound is enough: round up */
            mp_real_add(hi, hi, r, x->size + 2);
            hi->err = 0;
            mp_real_exp_bits(u, hi, 30);
            mp_real_zero(res);
            _mp_real_elem_add_mag(res, u);
            mp_real_clear(hi);
            mp_real_clear(r);
            mp_real_clear(u);
            return;
        }

        prec = FLINT_MIN(prec, acc + 8);
    }

    _ex_mid(res, &mid, prec);

    if (xerr != 0)
    {
        /* exp over the ball is exp(m) [1 +- u], u = e^rho - 1 <=
           rho (1 + rho) (rho < 1/4): the radius grows by
           |res| rho (1 + rho), in double arithmetic rounded up */
        double rd = (double) xerr;
        double rho, v;

        /* rho = rd B^xanc < 1/4, so xanc <= -1; clamped at 2^-896
           (normal) from below */
        FLINT_ASSERT(xanc <= -1);
        rho = _mp_real_elem_scale_up(rd, xanc);
        v = rd * (1.0 + rho) * _mp_real_elem_mag_hi(res) * (1.0 + 0x1p-48);
        _mp_real_elem_add_rad_d(res, v, xanc + res->exp - 1);
    }
}
