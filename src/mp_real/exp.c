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
_ex_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n, int notab)
{
    slong z = _mp_real_elem_lzb(v, n);

    if (z == WORD_MAX)
    {
        flint_mpn_zero(y, n);
        y[n] = 1;
        *err = 0;
    }
    else if (notab)
        _mp_real_exp_notab(y, err, v, n);
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
_ex_reduce(mp_real_t res, const mp_real_t m, slong n, int notab)
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
    _ex_kernel(res->d, &err, X + 1, n, notab);
    _mp_real_elem_finish(res, n, err + 4, 0);
    mp_real_mul_2exp_si(res, res, k);

    TMP_END;
}

/* exp(m), m exact, nonzero, |m| < 2^(FLINT_BITS - 5), to about prec
   bits */
static void
_ex_mid(mp_real_t res, const mp_real_t m, slong prec, int notab)
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
        _ex_kernel(res->d, &err, v, n, notab);
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

    _ex_reduce(res, m, n, notab);
}

/* ==== exp for huge precision without tables ================================

   Three ways to exp(m) for exact m != 0, |m| < 2^(FLINT_BITS - 5), none
   of which builds the tables of the bitwise or diophantine kernels (one
   evaluation at millions of bits does not pay for them):

   NOTAB_LOG2: the evaluation of mp_real_exp_bits (its short series for
   small |m|, else the reduction m = k log 2 + t with t in [0, log 2))
   with _mp_real_exp_notab as the kernel; log 2 at the working precision
   comes from (and fills) the constant cache.

   NOTAB_SQUARING: exp(|m|) = exp(|m| 2^-s)^(2^s) with |m| 2^-s < 1/2,
   _mp_real_exp_notab on the reduced argument and s ball squarings, one
   reciprocal for m < 0; no log 2.  The relative error doubles per
   squaring: s + 16 guard bits.  arb uses the same scheme above 10^6
   bits; it pays when log 2 is not already cached.

   AGM: one Newton-Taylor step on the AGM logarithm.  y0 = exp(m) to
   about P/(N + 1) limbs by NOTAB_SQUARING, taken as exact; with
   s = y0 2^j >= 2^(p/2 + 16) (j >= 0, p = FLINT_BITS P the working
   precision),

       log s = pi / (2 agm(1, 4/s)) + eps,  |eps| <= 64 (log s + 8) / s^2

   (Borwein & Borwein, Pi and the AGM, Thm. 7.2), so log y0 = log s - j
   log 2 (no log 2 at all when y0 is already large enough, as for the
   exp(C) of the partition function), d = m - log y0 with |d| < B^-(P /
   (N + 1)), and exp(m) = y0 exp(d), expm1(d) by the N-term Taylor
   series of the Newton steps of -log and atan (newton.c), coefficients
   N!/k! over the common denominator N!.  If y0 turns out too
   inaccurate for the series (not expected), the squaring result at full
   precision is used. */

#define EX_METHOD_AUTO 0
#define EX_METHOD_NOTAB_LOG2 1
#define EX_METHOD_SQUARING 2
#define EX_METHOD_AGM 3

/* |m| < 2^emid, |m| >= 2^(emid - 1) */
static slong
_ex_bitexp(const mp_real_t m)
{
    return FLINT_BITS * (m->exp - 1) + FLINT_BIT_COUNT(m->d[m->size - 1]);
}

static void
_ex_mid_squaring(mp_real_t res, const mp_real_t m, slong prec)
{
    slong emid = _ex_bitexp(m), s, nk, i;
    mp_real_t t;
    nn_ptr v;
    ulong er, e2;
    TMP_INIT;

    s = FLINT_MAX(0, emid + 1);
    /* the notab bound (128 ulps), the truncation of the argument, s
       doublings of the relative error and s + 1 roundings */
    nk = (prec + s + FLINT_BIT_COUNT(s + 1) + 24 + FLINT_BITS - 1) / FLINT_BITS;

    mp_real_init(t);
    TMP_START;
    v = TMP_ALLOC(nk * sizeof(ulong));

    /* v = |m| 2^-s in [0, 1/2), truncated (er ulps) */
    mp_real_mul_2exp_si(t, m, -s);
    t->negative = 0;
    _mp_real_get_fixed(v, &er, t, nk);

    mp_real_fit_length(res, nk + 1);
    _mp_real_exp_notab(res->d, &e2, v, nk);
    /* the truncation of v moves exp(v) < 2 by twice as many ulps */
    _mp_real_elem_finish(res, nk, e2 + 2 * er, 0);

    for (i = 0; i < s; i++)
        mp_real_mul(res, res, res, nk + 1);

    if (m->negative)
    {
        mp_real_set_ui(t, 1);
        mp_real_div(res, t, res, nk + 1);
    }

    TMP_END;
    mp_real_clear(t);
}

#if FLINT_BITS == 64
#define EX_AGM_MAX_N 16         /* 16! < 2^64 */
#else
#define EX_AGM_MAX_N 12         /* 12! < 2^32 */
#endif

static slong
_ex_agm_default_N(slong P)
{
    return FLINT_MIN((P <= 20000) ? 8 : 12, EX_AGM_MAX_N);
}

static void
_ex_mid_agm(mp_real_t res, const mp_real_t m, slong prec)
{
    slong P, N, w0, lg, j, T, uexp, k, emid = _ex_bitexp(m);
    mp_real_t y0, b, a, L, u, d;
    ulong c[EX_AGM_MAX_N], den;

    /* d = m - log y0 is computed as m - (log s - j log 2), terms of
       magnitude up to about |m| + p: their absolute errors need
       log2(|m| + p) guard bits on top of prec */
    P = (prec + FLINT_BITS - 1) / FLINT_BITS + 1;
    P += (FLINT_MAX(emid, 0) + FLINT_BIT_COUNT((ulong) (FLINT_BITS * P)) + 16) / FLINT_BITS + 1;

    N = _ex_agm_default_N(P);
    w0 = (P + N) / (N + 1) + 1;

    /* nothing to gain from the step at small precision */
    if (w0 + 2 >= P)
    {
        _ex_mid_squaring(res, m, prec);
        return;
    }

    mp_real_init(y0); mp_real_init(b); mp_real_init(a);
    mp_real_init(L); mp_real_init(u); mp_real_init(d);

    /* y0 ~ exp(m), exact */
    _ex_mid_squaring(y0, m, FLINT_BITS * w0);
    y0->err = 0;
    _mp_real_norm(y0);
    lg = _ex_bitexp(y0) - 1;                    /* y0 >= 2^lg */

    /* s = y0 2^j >= 2^T */
    T = FLINT_BITS * P / 2 + 16;
    j = FLINT_MAX(0, T - lg);

    /* L = pi / (2 agm(1, 4/s)) - j log 2 = log y0 + eps */
    mp_real_mul_2exp_si(u, y0, j);
    mp_real_set_ui(b, 4);
    mp_real_div(b, b, u, P);
    mp_real_set_ui(a, 1);
    mp_real_agm(b, a, b, P);            /* takes over b */
    mp_real_clear(a);
    mp_real_const_pi4(u, P, 1);
    mp_real_div(L, u, b, P);
    mp_real_clear(b);
    mp_real_mul_2exp_si(L, L, 1);
    /* 64 (log s + 8) / s^2 with log s < lg + j + 1 */
    mp_real_add_error_2exp_si(L, 6 + FLINT_BIT_COUNT((ulong) (lg + j + 9)) - 2 * (lg + j));
    if (j > 0)
    {
        mp_real_const_log2(u, P, 1);
        mp_real_mul_ui(u, u, (ulong) j, P);
        mp_real_sub(L, L, u, P);
    }

    /* d = m - log y0 (the full-length temporaries are freed as soon as
       they are dead: this is the peak of a huge-precision exp) */
    mp_real_sub(d, m, L, P);
    mp_real_clear(L);
    uexp = mp_real_abs_bound_lt_2exp_si(d);

    if ((N + 1) * (-uexp) < FLINT_BITS * P + 8)
    {
        /* y0 too inaccurate for the series (not expected) */
        _ex_mid_squaring(res, m, prec);
    }
    else
    {
        /* expm1(d) = d sum_{k=1}^N d^(k-1)/k! + tail, the tail below
           2 |d|^(N+1) / (N+1)! < 2^((N+1) uexp + 1) */
        den = 1;
        for (k = 2; k <= N; k++)
            den *= (ulong) k;
        c[N - 1] = 1;
        for (k = N - 1; k >= 1; k--)
            c[k - 1] = c[k] * (ulong) (k + 1);      /* N!/k! */
        _mp_real_newton_series(u, d, uexp, P, c, N, 0, 0);
        mp_real_div_ui(u, u, den, P);
        mp_real_add_error_2exp_si(u, (N + 1) * uexp + 1);
        mp_real_clear(d);
        mp_real_init(d);
        mp_real_set_ui(d, 1);           /* d no longer needed: 1 */
        mp_real_add(u, u, d, P);
        mp_real_mul(res, y0, u, P);
    }

    mp_real_clear(y0); mp_real_clear(u); mp_real_clear(d);
}

static void
_ex_mid_method(mp_real_t res, const mp_real_t m, slong prec, int method)
{
    if (method == EX_METHOD_NOTAB_LOG2)
        _ex_mid(res, m, prec, 1);
    else if (method == EX_METHOD_SQUARING)
        _ex_mid_squaring(res, m, prec);
    else if (method == EX_METHOD_AGM)
        _ex_mid_agm(res, m, prec);
    else
        _ex_mid(res, m, prec, 0);
}

static void
_ex_ball(mp_real_t res, const mp_real_t x, slong prec, int method)
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
        flint_throw(FLINT_ERROR, "mp_real_exp: |x| >= 2^%wd, the result "
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

    _ex_mid_method(res, &mid, prec, method);

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

void
mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    _ex_ball(res, x, prec, EX_METHOD_AUTO);
}

void
mp_real_exp_notab_log2(mp_real_t res, const mp_real_t x, slong n)
{
    _ex_ball(res, x, FLINT_BITS * FLINT_MAX(n, 1), EX_METHOD_NOTAB_LOG2);
}

void
mp_real_exp_notab_squaring(mp_real_t res, const mp_real_t x, slong n)
{
    _ex_ball(res, x, FLINT_BITS * FLINT_MAX(n, 1), EX_METHOD_SQUARING);
}

void
mp_real_exp_agm(mp_real_t res, const mp_real_t x, slong n)
{
    _ex_ball(res, x, FLINT_BITS * FLINT_MAX(n, 1), EX_METHOD_AGM);
}
