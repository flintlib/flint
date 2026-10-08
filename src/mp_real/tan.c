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

/* tan of a ball to a relative accuracy of about 2^-prec.

   The midpoint m is evaluated exactly as given, by the fixed-point tan
   kernel on |m| for |m| < 1 (at enough limbs that the result keeps its
   relative accuracy; small arguments take a short Taylor polynomial
   instead) and otherwise after the reduction of mp_real_sin_cos_bits,
   |m| = a pi/2 + s v with v in [0, pi/4]:

       tan m = s tan v (a even),   -s / tan v (a odd).

   The reduction is absolute (3 ulps), and near a zero or a pole of the
   tangent v is small, so an absolute error is a relative one of tan m:
   the reduction is repeated with as many more limbs as v has leading
   zero limbs, until v is accurate to the working precision relative to
   itself.  An exact m is never a multiple of pi/2, and by the
   irrationality measure of pi (below 7.2) a b-bit m lies at least about
   2^(-7.2 b) away from one: the repetition is capped at 8 b + 1024 more
   bits, beyond which the status is 0.

   The radius r of x does not enter as for the sine and cosine, whose
   derivatives are bounded: with T enclosing tan m and U = r (1 + |T|),
   which is at least r / |cos m|, U < 1 means that |cos| (1-Lipschitz)
   stays above |cos m| - r > 0 over the ball, which then contains no pole,
   and

       |tan y - tan m| <= r max sec^2 <= r / (|cos m| - r)^2
                        <= r (1 + T^2) / (1 - U)^2,

   evaluated in two-limb balls (_mp_real_tan_add_arg_rad).  Otherwise the
   status is 0: the ball may contain a pole (the test is conservative
   by a factor below sqrt(2), where the output would have almost no
   accuracy anyway). */

/* kernel choice by working precision (limbs), as for sin and cos */
#define TAN_BITWISE_MAX 600
#define TAN_DIOPHANTINE_MAX 65536

/* the leading zero bits from which sin and 1 - cos of the reduced
   argument (_mp_real_sin_cos_reduced) and a division beat the tan
   kernel, where the tan series (below) does not apply (measured: the
   pair costs 0.6-0.9 of the kernel there, the division about 0.1) */
static slong
_tan_reduced_min_z(slong n)
{
    if (n <= 32)
        return 32;
    if (n <= 64)
        return 40;
    if (n <= 96)
        return 48;
    if (n <= 192)
        return 64;
    if (n <= 256)
        return 96;
    if (n <= TAN_BITWISE_MAX)
        return 160;
    if (n <= 1024)
        return 128;
    if (n <= 4096)
        return 160;
    if (n <= 8192)
        return 256;
    if (n <= TAN_DIOPHANTINE_MAX)
        return 4096;
    return WORD_MAX;    /* the notab kernel adapts by itself */
}

/* (y, n + 1) with a units limb from a ball enclosing a value in [0, 2)
   (clamped at 0 from below); returns the bound in ulps of B^-n,
   UWORD_MAX when none fits a word */
ulong
_mp_real_ball_get_fixed1(nn_ptr y, const mp_real_t x, slong n)
{
    ulong e;

    if (x->size != 0 && !x->negative && x->exp >= 1)
    {
        slong k = x->exp - x->size + n;
        int dropped;

        if (x->exp > 1)
            return UWORD_MAX;
        dropped = _mp_real_elem_copy(y, n + 1, n, x);
        if (x->err == 0)
            e = 0;
        else if (k < 0)
            e = 1;      /* err < B units of B^(k - n) < one ulp */
        else if (k == 0)
            e = x->err;
        else
            return UWORD_MAX;
        return (e == UWORD_MAX) ? e : e + dropped;
    }

    _mp_real_get_fixed(y, &e, x, n);
    y[n] = 0;
    return e;
}

/* (y, n + 1) = s / c for s in [0, 1), c in (1/2, 1] given as fixed-point
   balls (s, n + 1) +- es, (c, n + 1) +- ec ulps of B^-n */
static ulong
_tan_div(nn_ptr y, nn_srcptr s, ulong es, nn_srcptr c, ulong ec, slong n)
{
    mp_real_t S, C;
    ulong e;

    mp_real_init(S);
    mp_real_init(C);
    mp_real_fit_length(S, n + 1);
    mp_real_fit_length(C, n + 1);
    flint_mpn_copyi(S->d, s, n + 1);
    flint_mpn_copyi(C->d, c, n + 1);
    _mp_real_elem_finish(S, n, FLINT_MAX(es, 1), 0);
    _mp_real_elem_finish(C, n, FLINT_MAX(ec, 1), 0);
    mp_real_div(S, S, C, n + 1);
    e = _mp_real_ball_get_fixed1(y, S, n);
    mp_real_clear(S);
    mp_real_clear(C);
    return e;
}

/* forced inline into its callers here (see sin_cos.c), exported
   through the wrapper below */
FLINT_FORCE_INLINE void
_tan_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n)
{
    slong z = _mp_real_elem_lzb(v, n);
    ulong e;

    if (z == WORD_MAX)
    {
        flint_mpn_zero(y, n + 1);
        e = 0;
    }
    else if (z >= MP_REAL_SERIES_TAN_RMIN
        && _mp_real_series_rs_tan_ok(n, (flint_bitcnt_t) z))
    {
        /* the tangent series: faster than every kernel (0.1-0.9 of the
           kernel's time where it applies) */
        _mp_real_series_rs_tan(y, v, n, (flint_bitcnt_t) z);
        y[n] = 0;
        e = 3;
    }
    else if (z >= 32 && z >= _tan_reduced_min_z(n))
    {
        /* sin v and g = 1 - cos v of the reduced argument, divided */
        nn_ptr s, c;
        ulong es;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(2 * (n + 1) * sizeof(ulong));
        c = s + n + 1;
        _mp_real_sin_cos_reduced(s, c, &es, v, n, (flint_bitcnt_t) z, 0);
        s[n] = 0;
        c[n] = mpn_neg(c, c, n) ? 0 : 1;
        e = _tan_div(y, s, es, c, es, n);
        TMP_END;
    }
    else if (n <= TAN_BITWISE_MAX)
        _mp_real_tan_bitwise_rs(y, &e, v, n, 0);
    else if (n <= TAN_DIOPHANTINE_MAX)
        _mp_real_tan_diophantine(y, &e, v, n);
    else
    {
        /* no untabulated tangent: sin and cos, divided */
        nn_ptr s, c;
        ulong es;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(2 * (n + 1) * sizeof(ulong));
        c = s + n + 1;
        _mp_real_sin_cos_notab(s, c, &es, v, n);
        e = _tan_div(y, s, es, c, es, n);
        TMP_END;
    }

    if (err != NULL)
        *err = e;
}

void
_mp_real_tan_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n)
{
    _tan_kernel(y, err, v, n);
}

/* guard bits for a kernel result accurate to about 2^-p (absolute):
   the kernel's bound (8r + 256 for the bitwise functions, r the tuned
   parameter; small for the others), a reduction error of a few ulps
   amplified by sec^2 <= 2, and two bits of slack */
static slong
_tan_guard(slong p)
{
    slong n = (p + 14 + FLINT_BITS - 1) / FLINT_BITS;

    if (n <= TAN_BITWISE_MAX)
        return FLINT_BIT_COUNT(8 * (ulong) _mp_real_trig_default_r_inline(n) + 264) + 3;
    if (n <= TAN_DIOPHANTINE_MAX)
        return 10;
    return 14;
}

slong
_mp_real_tan_limbs(slong p)
{
    slong n = (p + _tan_guard(p) + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* res = the lower resp. upper endpoint of a strictly positive ball x,
   exactly (the radius is a count of ulps of x's bottom limb) */
static void
_ball_endpoint(mp_real_t res, const mp_real_t x, int upper)
{
    slong n = x->size;
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC((n + 1) * sizeof(ulong));
    flint_mpn_copyi(t, x->d, n);
    t[n] = 0;
    if (upper)
        t[n] = mpn_add_1(t, t, n, x->err);
    else
        mpn_sub_1(t, t, n, x->err);
    _mp_real_set_mpn_2exp(res, t, n + 1, FLINT_BITS * (x->exp - n));
    TMP_END;
}

int
_mp_real_ball_inv_pos(mp_real_t res, const mp_real_t x, slong n)
{
    mp_real_t one;

    if (x->size == 0 || x->negative || !(_mp_real_elem_mag_lo(x) > 0.0))
        return 0;

    mp_real_init(one);
    mp_real_set_ui(one, 1);

    if (mp_real_rel_radius_lt_2exp_si(x) <= -32)
        mp_real_div(res, one, x, n);
    else
    {
        /* [1/hi, 1/lo] from the exact endpoints, as the ball of their
           midpoint and half-difference */
        mp_real_t lo, hi;

        mp_real_init(lo);
        mp_real_init(hi);
        _ball_endpoint(lo, x, 0);
        _ball_endpoint(hi, x, 1);
        mp_real_div(lo, one, lo, n);
        mp_real_div(hi, one, hi, n);
        mp_real_add(res, lo, hi, n);
        mp_real_sub(lo, lo, hi, n);
        mp_real_mul_2exp_si(res, res, -1);
        mp_real_mul_2exp_si(lo, lo, -1);
        _mp_real_elem_add_mag(res, lo);
        mp_real_clear(lo);
        mp_real_clear(hi);
    }

    mp_real_clear(one);
    return 1;
}

int
_mp_real_inv_mpn(mp_real_t res, nn_srcptr Y, slong L, slong ey, ulong err,
    slong rx, slong wn, int neg)
{
    slong an, k, bt, ki[4], i;
    double S;
    nn_ptr q;
    TMP_INIT;

    while (L > 0 && Y[L - 1] == 0)
        L--;
    if (L == 0)
        return 0;
    bt = FLINT_BIT_COUNT(Y[L - 1]);

    /* the relative errors, 2^ki bounding each: y (err / Y, with
       Y >= 2^(bt - 1) B^(L - 1)), the further factor 2^rx, the
       truncation of Y to its top an limbs and the Newton inverse (4 B^-wn
       relative, the inverse of a in [1/B, 1) exceeding 1) */
    an = FLINT_MIN(L, wn + 2);
    ki[0] = (err != 0) ? (slong) FLINT_BIT_COUNT(err) - (bt - 1) - FLINT_BITS * (L - 1) : WORD_MIN;
    ki[1] = rx;
    ki[2] = (an < L) ? -(bt - 1) - FLINT_BITS * (an - 1) : WORD_MIN;
    ki[3] = 2 - FLINT_BITS * wn;
    k = ki[3];
    for (i = 0; i < 3; i++)
        k = FLINT_MAX(k, ki[i]);

    /* the reciprocal of a ball with relative radius below 2^-12 */
    if (k > -12)
    {
        mp_real_t b;
        int ok;

        mp_real_init(b);
        _mp_real_set_mpn_2exp(b, Y, L, FLINT_BITS * ey);
        if (err != 0)
            _mp_real_elem_add_rad(b, err, ey);
        if (rx != WORD_MIN)
            mp_real_add_error_2exp_si(b, rx + 1 + mp_real_abs_bound_lt_2exp_si(b));
        ok = _mp_real_ball_inv_pos(b, b, wn);
        if (ok)
        {
            if (neg)
                mp_real_neg(b, b);
            mp_real_swap(res, b);
        }
        mp_real_clear(b);
        return ok;
    }

    TMP_START;
    q = TMP_ALLOC((wn + 2) * sizeof(ulong));

    /* q = 1/a within 4 B^-wn for a = (Y top an limbs) B^-an, and
       1/y = (1/a) B^(-ey - L) */
    _mp_real_inv_newton(q, Y + (L - an), an, wn);
    _mp_real_set_mpn_2exp(res, q, wn + 2, -FLINT_BITS * (wn + ey + L));
    if (neg)
        res->negative = 1;

    /* relatively, 1/y_true deviates from res by at most
       (1 + e3)(1 + e2) / ((1 - e0)(1 - e1)) - 1 <= 1.001 (e0 + ... + e3)
       = S 2^k for all ei <= 2^-12, and |1/y_true| <= |res| (1 + 2^-10):
       an absolute bound 1.002 S 2^k |res|, |res| < (top + 1) B^(exp - 1) */
    S = 0.0;
    for (i = 0; i < 4; i++)
        if (ki[i] != WORD_MIN && ki[i] - k > -60)
            S += d_mul_2exp_inrange(1.0, (int) (ki[i] - k));
        else if (ki[i] != WORD_MIN)
            S += 0x1p-60;
    S *= 1.0021;
    {
        /* S 2^k (top + 1) B^(exp - 1), as v B^(exp - 1 + floor(k / B))
           with v = S (top + 1) 2^(k mod B), in [2^-64, 2^66) */
        slong kq = k >> MP_REAL_LGB;
        int kr = (int) (k - kq * FLINT_BITS);
        double v = S * ((double) res->d[res->size - 1] + 1.0) * d_mul_2exp_inrange(1.0, kr - FLINT_BITS);
        _mp_real_elem_add_rad_d(res, v * (1.0 + 0x1p-50), res->exp - 1 + kq + 1);
    }
    TMP_END;
    return 1;
}

int
_mp_real_tan_add_arg_rad(mp_real_t res, const mp_real_t R)
{
    mp_real_t T, Y, U, D, E, one;
    int ok = 1;

    if (R->size == 0 && R->err == 0)
        return 1;

    mp_real_init(T);
    mp_real_init(Y);
    mp_real_init(U);
    mp_real_init(D);
    mp_real_init(E);
    mp_real_init(one);
    mp_real_set_ui(one, 1);

    /* T encloses |tan m| at two limbs: |res|, truncated, with its radius
       (||t| - |mid|| <= |t - mid|) */
    if (res->size == 0)
        mp_real_set(T, res);
    else
    {
        _mp_real_elem_set_trunc(T, res, 2);
        if (res->err != 0)
            _mp_real_elem_add_rad(T, res->err, res->exp - res->size);
    }
    T->negative = 0;

    /* U = R (1 + |T|) >= r / |cos m|; D = 1 - U must be positive */
    mp_real_add(U, T, one, 2);
    mp_real_mul(U, U, R, 2);
    mp_real_sub(D, one, U, 2);

    if (D->size == 0 || D->negative || !(_mp_real_elem_mag_lo(D) > 0.0))
    {
        ok = 0;
    }
    else
    {
        /* E = R (1 + T^2) / D^2, bounded through the lower endpoint of D
           (exact, so that the division takes any width of D) */
        mp_real_mul(Y, T, T, 2);
        mp_real_add(Y, Y, one, 2);
        mp_real_mul(E, R, Y, 2);
        _ball_endpoint(U, D, 0);
        mp_real_mul(D, U, U, 2);
        mp_real_div(E, E, D, 2);
        _mp_real_elem_add_mag(res, E);
    }

    mp_real_clear(T);
    mp_real_clear(Y);
    mp_real_clear(U);
    mp_real_clear(D);
    mp_real_clear(E);
    mp_real_clear(one);
    return ok;
}

/* res = tan x for an exact x with |x| < 2^e, 2 (-e) >= prec + 3: x
   within 2^(3e - 1) (tan y = y + c_1 y^3 + ..., c_1 = 1/3 and
   c_(k+1) / c_k <= 0.41, below 0.372 |y|^3 for |y| <= 1/2).  Beyond, the
   kernel's series for small arguments beats a polynomial in ball
   arithmetic (measured). */
static int
_tan_small(mp_real_t res, const mp_real_t x, slong prec)
{
    slong e = mp_real_abs_bound_lt_2exp_si(x);

    if (e > -1 || -2 * e < prec + 3)
        return 0;

    if (res != x)
        mp_real_set(res, x);
    mp_real_add_error_2exp_si(res, 3 * e - 1);
    return 1;
}

/* res = the ball (y, n + 1) B^-n +- err ulps (err >= 1), sign neg */
static void
_tan_set_fixed(mp_real_t res, nn_srcptr y, slong n, ulong err, int neg)
{
    mp_real_fit_length(res, n + 1);
    flint_mpn_copyi(res->d, y, n + 1);
    _mp_real_elem_finish(res, n, FLINT_MAX(err, 1), neg);
}

/* |x| < 1, exact, nonzero, not small: the kernel on the fixed-point |x|
   at n fraction limbs (x truncated there, one more ulp, times
   sec^2 x <= 3.43) */
static void
_tan_unit(mp_real_t res, const mp_real_t x, slong n)
{
    nn_ptr buf, y;
    nn_srcptr v;
    ulong err;
    int trunc, xneg = x->negative;
    TMP_INIT;

    TMP_START;
    buf = TMP_ALLOC(n * sizeof(ulong));
    y = TMP_ALLOC((n + 1) * sizeof(ulong));
    v = _mp_real_elem_frame(buf, n, x, res->d, NULL, &trunc);
    _tan_kernel(y, &err, v, n);
    _tan_set_fixed(res, y, n, err + 4 * trunc, xneg);
    TMP_END;
}

/* |x| >= 1, exact: reduce mod pi/2 at the limbs that make v (z leading
   zero bits) relatively accurate, tan_limbs(prec + z), repeating the
   reduction when z calls for more than the first pass had, and evaluate
   into res (which may alias x: x is read by the reductions only); returns
   0 if a pole cannot be excluded.  noretry_z: the leading zero bits of v
   from which no repetition pays (the argument's radius dominates) */
static int
_tan_reduce(mp_real_t res, const mp_real_t x, slong prec, slong noretry_z)
{
    slong N0, N, Nn, z, wn;
    slong cap = 8 * x->size + 1024 / FLINT_BITS + 1;
    nn_ptr W, y;
    ulong err;
    int code, a, s, ok = 1, neg;
    TMP_INIT;

    TMP_START;
    N = N0 = _mp_real_tan_limbs(prec + 2);
    for (;;)
    {
        W = TMP_ALLOC(N * sizeof(ulong));
        code = _mp_real_trig_reduce(W, x, N);
        z = _mp_real_elem_lzb(W, N);

        if (z == WORD_MAX)
        {
            if (FLINT_BITS * N >= noretry_z)
                break;
            Nn = 2 * N;
        }
        else
        {
            /* v within 3 ulps: relatively accurate to tan_limbs' target
               when FLINT_BITS N covers prec + z and the guard */
            if (z >= noretry_z)
                break;
            Nn = _mp_real_tan_limbs(prec + z);
            if (Nn <= N)
                break;
        }

        if (Nn - N0 > cap)
        {
            TMP_END;
            return 0;
        }
        N = Nn;
    }

    a = code & 3;
    s = (code & 4) ? -1 : 1;

    /* tan m = s tan v (a even), -s / tan v (a odd), negated for m < 0 */
    neg = ((a & 1) ? (s > 0) : (s < 0)) ^ x->negative;

    /* tan v within the kernel's bound plus the reduction's 3 ulps times
       sec^2 v <= 2 (1 + 2/B) */
    y = TMP_ALLOC((N + 1) * sizeof(ulong));
    _tan_kernel(y, &err, W, N);
    err += 7;

    if (!(a & 1))
        _tan_set_fixed(res, y, N, err, neg);
    else
    {
        wn = mp_real_prec_bits(prec);
        ok = _mp_real_inv_mpn(res, y, N + 1, -N, err, WORD_MIN, wn, neg);
    }

    TMP_END;
    return ok;
}

int
mp_real_tan_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    slong e, noretry_z = WORD_MAX;
    ulong xerr = x->err;
    slong xanc = x->exp - x->size;
    int ok = 1;

    prec = FLINT_MAX(prec, 2);

    /* x = 0 exactly */
    if (x->size == 0 && x->err == 0)
    {
        mp_real_zero(res);
        return 1;
    }

    /* |x| < 2^e over the whole ball */
    e = mp_real_abs_bound_lt_2exp_si(x);

    /* a ball around zero, or one too wide for a relative accuracy of two
       bits, within (-1, 1): |tan y| <= tan(1) |y| < 2^(e + 1) */
    if (e <= 0 && (x->size == 0 || (x->err != 0
            && mp_real_rel_radius_lt_2exp_si(x) > -2)))
    {
        mp_real_zero(res);
        mp_real_add_error_2exp_si(res, e + 1);
        return 1;
    }

    /* ridiculously large arguments */
    if (FLINT_BITS * (x->exp - 1) >= FLINT_MAX(65536, 4 * prec))
    {
        mp_real_zero(res);
        return 0;
    }

    /* the midpoint, exactly */
    mid = *x;
    mid.err = 0;

    /* an inexact x: the output's relative error is at least that of x
       for |x| < 1 (tan y / y >= 1, sec^2 y >= 1) and at least 2r, r the
       absolute radius, otherwise ((1 + t^2) / |t| >= 2) */
    if (xerr != 0 && x->size != 0)
    {
        slong rel = mp_real_rel_radius_lt_2exp_si(x), acc;

        if (x->exp <= 0)
            acc = -rel;
        else
            acc = -(rel + FLINT_BITS * (x->exp - 1) + FLINT_BIT_COUNT(x->d[x->size - 1]));

        prec = FLINT_MIN(prec, FLINT_MAX(acc + 8, 10));

        /* near a pole or zero, v below the absolute radius makes no
           repeated reduction worthwhile */
        noretry_z = FLINT_MAX(-(FLINT_BITS * xanc + (slong) FLINT_BIT_COUNT(xerr)) + 8, 0);
    }

    /* (the output may alias x: from here on, only mid, xerr and xanc are
       read, and each path reads mid before writing res) */
    if (mid.size == 0)
    {
        mp_real_zero(res);
    }
    else if (mid.exp <= 0)
    {
        /* |x| < 2^emid, emid <= 0 */
        slong emid = FLINT_BITS * (mid.exp - 1) + FLINT_BIT_COUNT(mid.d[mid.size - 1]);
        slong z = -emid;

        if (!_tan_small(res, &mid, prec))
            _tan_unit(res, &mid, _mp_real_tan_limbs(prec + z));
    }
    else
    {
        ok = _tan_reduce(res, &mid, prec, noretry_z);
    }

    /* the radius of x */
    if (ok && xerr != 0)
    {
        mp_real_t R;
        mp_real_init(R);
        mp_real_set_ui(R, xerr);
        mp_real_mul_2exp_si(R, R, FLINT_BITS * xanc);
        ok = _mp_real_tan_add_arg_rad(res, R);
        mp_real_clear(R);
    }

    if (!ok)
        mp_real_zero(res);

    return ok;
}
