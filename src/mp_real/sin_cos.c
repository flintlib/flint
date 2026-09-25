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

/* sin and cos of a ball to a relative accuracy of about 2^-prec.

   The midpoint is evaluated exactly as given; the radius r of x enters
   the outputs at the end (|sin'|, |cos'| <= 1), and prec is first
   lowered to what r allows.  For |x| < 1 the argument goes to the
   fixed-point kernel as it is, at enough limbs that sin keeps its
   relative accuracy (small arguments take a short Taylor polynomial
   instead).  For |x| >= 1 it is reduced mod pi/2: with q a lower
   approximation of x 2/pi (from a cached 2/pi, never above the exact
   quotient) and P <= pi/2 the cached pi/2 at m + n + 1 fraction limbs,
   t = x - q P lies in [0, pi/2 (1 + 2/B)).  For t > pi/4 the angle is
   folded, v = P - t (slightly negative only when q came out one low).
   So x = a pi/2 + s v, v in [0, pi/4], and

       sin x = s sin v, cos v, -s sin v, -cos v   for a = 0, 1, 2, 3,
       cos x = cos v, -s sin v, -cos v, s sin v,

   where the kernel gives sin v and cos v at n fraction limbs.  The
   reduction errs by q (pi/2 - P) + truncation < 3 ulps of B^-n. */

/* kernel choice by working precision (limbs) */
#define SC_BITWISE_MAX 600
#define SC_DIOPHANTINE_MAX 65536

/* the series from z leading zero bits on: sin and 1 - cos by the
   reduced series of series_tapered.c / series_rs.c below 32 zero bits
   (and below 48 at 11-16 limbs), _mp_real_sin_cos_reduced beyond
   (measured end to end, kernel / series time ratios: 1.1-1.4 at 2-8
   limbs from z = 18, 1.1-1.3 at 10-20 limbs from z = 24-28, 1.1-1.6 at
   24-128 limbs from z = 32, 1.1 at 192-384 limbs from z = 40-64) */
static slong
_sc_series_min_z(slong n)
{
    if (n <= 1)
        return WORD_MAX;
    if (n <= 8)
        return 18;
    if (n <= 12)
        return 24;
    if (n <= 20)
        return 28;
    if (n <= 128)
        return 32;
    if (n <= 192)
        return 40;
    if (n <= 256)
        return 48;
    if (n <= 384)
        return 64;
    if (n <= SC_BITWISE_MAX)
        return 96;
    if (n <= 2048)
        return 160;
    if (n <= 8192)
        return 256;
    if (n <= SC_DIOPHANTINE_MAX)
        return 4096;
    return WORD_MAX;    /* the notab kernel adapts by itself */
}

/* (ys, n + 1), (yc, n + 1) = sin v, cos v for v in [0, 1) at n
   fraction limbs */
static void
_sc_kernel(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v, slong n)
{
    slong top = n - 1, z;

    while (top >= 0 && v[top] == 0)
        top--;
    z = (top < 0) ? WORD_MAX : FLINT_BITS * (n - 1 - top) + flint_clz(v[top]);

    if (z >= _sc_series_min_z(n)
        && !(z < 32 || (n >= 11 && n <= 16 && z < 48)))
    {
        if (z == WORD_MAX)
        {
            flint_mpn_zero(ys, n + 1);
            flint_mpn_zero(yc, n);
            yc[n] = 1;
            *err = 0;
        }
        else
        {
            /* sin v and g = 1 - cos v (below 1/2); z >= 32 */
            _mp_real_sin_cos_reduced(ys, yc, err, v, n, (flint_bitcnt_t) z, 0);
            ys[n] = 0;
            yc[n] = mpn_neg(yc, yc, n) ? 0 : 1;
        }
    }
    else if (z >= _sc_series_min_z(n))
    {
        /* sin v and g = 1 - cos v by the tapered series */
        _mp_real_series_sin_cos(ys, yc, v, n, (flint_bitcnt_t) z);
        ys[n] = 0;
        yc[n] = mpn_neg(yc, yc, n) ? 0 : 1;
        *err = 6;
    }
    else if (n <= SC_BITWISE_MAX)
        _mp_real_sin_cos_bitwise_rs(ys, yc, err, v, n, 0);
    else if (n <= SC_DIOPHANTINE_MAX)
        _mp_real_sin_cos_diophantine(ys, yc, err, v, n);
    else
        _mp_real_sin_cos_notab(ys, yc, err, v, n);
}

/* guard bits for a kernel result accurate to about 2^-p (absolute):
   the kernel's bound (6r + 128 for the bitwise functions, r the tuned
   parameter; at most 128 otherwise), the reduction's 3 ulps and two
   bits of slack */
static slong
_sc_guard(slong p)
{
    slong n = (p + 12 + FLINT_BITS - 1) / FLINT_BITS;

    if (n <= SC_BITWISE_MAX)
        return FLINT_BIT_COUNT(6 * (ulong) _mp_real_trig_bitwise_rs_default_r(n) + 131) + 3;
    return 10;
}

/* the kernel limbs for a target of about 2^-p (absolute) */
static slong
_sc_limbs(slong p)
{
    slong n = (p + _sc_guard(p) + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* sin x and cos x from v in [0, 1) at n fraction limbs, for
   x = a pi/2 + s v (s = +-1) up to eps ulps of B^-n, with the sine
   negated for a negative x.  v must not alias the outputs' limbs. */
static void
_sc_eval(mp_real_t rs, mp_real_t rc, nn_srcptr v, slong n, ulong eps,
    int a, int s, int xneg)
{
    /* sin v goes to rs for even a, to rc for odd a */
    mp_real_ptr os = (a & 1) ? rc : rs, oc = (a & 1) ? rs : rc;
    int sneg = (s < 0), negs, negc;
    nn_ptr ys, yc, tmp = NULL;
    ulong err;
    TMP_INIT;

    TMP_START;
    if (os == NULL || oc == NULL)
        tmp = TMP_ALLOC((n + 1) * sizeof(ulong));

    if (os != NULL)
    {
        mp_real_fit_length(os, n + 1);
        ys = os->d;
    }
    else
        ys = tmp;

    if (oc != NULL)
    {
        mp_real_fit_length(oc, n + 1);
        yc = oc->d;
    }
    else
        yc = tmp;

    _sc_kernel(ys, yc, &err, v, n);
    err += eps;

    /* sin x = s sin v, cos v, -s sin v, -cos v;
       cos x = cos v, -s sin v, -cos v, s sin v   (a = 0, 1, 2, 3) */
    if (a & 1)
    {
        negs = (a == 3);
        negc = sneg ^ (a == 1);
    }
    else
    {
        negs = sneg ^ (a == 2);
        negc = (a == 2);
    }

    if (rs != NULL)
        _mp_real_elem_finish(rs, n, err, negs ^ xneg);
    if (rc != NULL)
        _mp_real_elem_finish(rc, n, err, negc);

    TMP_END;
}

/* sin, cos of a small exact x, |x| < 2^e: sin x = x, cos x = 1 up to
   x^3/6 < 2^(3e-2) and x^2/2 = 2^(2e-1) */
static void
_sc_small1(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong e)
{
    if (rs != NULL)
    {
        mp_real_set(rs, x);
        rs->err = 0;
        mp_real_add_error_2exp_si(rs, 3 * e - 2);
    }
    if (rc != NULL)
    {
        mp_real_set_ui(rc, 1);
        mp_real_add_error_2exp_si(rc, 2 * e - 1);
    }
}

/* sin, cos of a small exact x, |x| < 2^e <= 1/2, as x - x^3/6 and
   1 - x^2/2 (remainders below x^5/120 < 2^(5e-6) and x^4/24 < 2^(4e-4))
   at w fraction limbs */
#define SC_SMALL2_FIXED_MAX 8

static void
_sc_small2(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong e, slong w,
    slong wp)
{
    if (w <= SC_SMALL2_FIXED_MAX)
    {
        /* fixed point: v = |x| truncated (1 ulp), x2 = floor(v^2) within
           2 ulps of x^2, x3 = floor(x2 v) within 2, x3/6 within 2, so
           sin within 3 ulps; x2/2 within 2, so cos within 2 */
        ulong v[SC_SMALL2_FIXED_MAX], t[2 * SC_SMALL2_FIXED_MAX];
        ulong x2[SC_SMALL2_FIXED_MAX];
        int xneg = x->negative;

        _mp_real_elem_copy(v, w, w, x);
        flint_mpn_sqr(t, v, w);
        flint_mpn_copyi(x2, t + w, w);

        if (rs != NULL)
        {
            flint_mpn_mul_n(t, x2, v, w);
            mpn_divrem_1(t + w, 0, t + w, w, 6);
            mp_real_fit_length(rs, w + 1);
            mpn_sub_n(rs->d, v, t + w, w);
            rs->d[w] = 0;
            _mp_real_elem_finish(rs, w, 3, xneg);
            mp_real_add_error_2exp_si(rs, 5 * e - 6);
        }

        if (rc != NULL)
        {
            mp_real_fit_length(rc, w + 1);
            mpn_rshift(t, x2, w, 1);
            rc->d[w] = 1 - mpn_neg(rc->d, t, w);
            _mp_real_elem_finish(rc, w, 2, 0);
            mp_real_add_error_2exp_si(rc, 4 * e - 4);
        }
    }
    else
    {
        mp_real_t x2, t;

        mp_real_init(x2);
        mp_real_init(t);
        mp_real_mul(x2, x, x, wp);

        if (rs != NULL)
        {
            mp_real_mul(t, x2, x, wp);
            mp_real_div_ui(t, t, 6, wp);
            mp_real_sub(t, x, t, wp);
            mp_real_add_error_2exp_si(t, 5 * e - 6);
            mp_real_swap(rs, t);
        }
        if (rc != NULL)
        {
            mp_real_mul_2exp_si(x2, x2, -1);
            mp_real_set_ui(t, 1);
            mp_real_sub(rc, t, x2, wp);
            mp_real_add_error_2exp_si(rc, 4 * e - 4);
        }

        mp_real_clear(x2);
        mp_real_clear(t);
    }
}

/* |x| < 1, exact, nonzero: the kernel on the fixed-point x at n
   fraction limbs (x truncated there, one more ulp) */
static void
_sc_unit(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong n)
{
    nn_ptr v;
    int trunc, xneg = x->negative;
    TMP_INIT;

    TMP_START;
    v = TMP_ALLOC(n * sizeof(ulong));
    trunc = _mp_real_elem_copy(v, n, n, x);
    _sc_eval(rs, rc, v, n, trunc, 0, 1, xneg);
    TMP_END;
}

#if FLINT_BITS == 64

/* the fraction limbs of pi/2 (low first), 2/pi = (R1, R0) / B^2
   truncated, and the top fraction limb of pi/4 */
#define SC_PI2_LIMBS 40
static const ulong _sc_pi2_frac[SC_PI2_LIMBS] = {
    UWORD(0xb0ec04e67d90d4c8), UWORD(0xe25ff40db31410c9), UWORD(0x9dc7a44c35a5dcd7),
    UWORD(0x3d1929c0944ac33b), UWORD(0x57eb5d19b61267ae), UWORD(0x672e1f0b4dc3c98f),
    UWORD(0x15d4e2aeba0c18fb), UWORD(0xd9f70a08b1b7de15), UWORD(0x50aa4357be3974c9),
    UWORD(0x5a662e1a08a0f467), UWORD(0x2ae51cb51555885b), UWORD(0x2ba44c3131f40a20),
    UWORD(0x732a92f9d52ad5ca), UWORD(0xbc5797ed2ab02e30), UWORD(0x6b8abbe0de98a593),
    UWORD(0x364f0745d80f451f), UWORD(0xc73cee58301d0c07), UWORD(0x6520bc8c5c6d9c77),
    UWORD(0xe2e8d811943042f8), UWORD(0xce186a9c95793009), UWORD(0x3daa520ee12d2cda),
    UWORD(0x38c5e6ac410aa577), UWORD(0x06caba47b9475b2c), UWORD(0xd22c7f51fa499ebf),
    UWORD(0x31b4906c38aba734), UWORD(0x8400f97142c77e0b), UWORD(0x9250cca3d9c8b67b),
    UWORD(0x5d3e4822f8963fcc), UWORD(0xdc70d7f6b5133f4b), UWORD(0x17feb96de80d6fdb),
    UWORD(0xe89885d34c6fdad6), UWORD(0xc90b6aecc4bcfd8d), UWORD(0x9fc26adadaa3848b),
    UWORD(0x605614dbe4be286e), UWORD(0xdf2a33679a748636), UWORD(0xa29410f31c6809bb),
    UWORD(0x04177d4c76273644), UWORD(0x52049c1114cf98e8), UWORD(0x898cc51701b839a2),
    UWORD(0x921fb54442d18469)
};
#define SC_2_DIV_PI_1 UWORD(0xa2f9836e4e441529)
#define SC_2_DIV_PI_0 UWORD(0xfc2757d1f534ddc0)
#define SC_PI4_TOP UWORD(0xc90fdaa22168c234)

/* X[0, N) -= q P[0, N), returning the borrow limb; P, q in registers
   for constant N */
FLINT_FORCE_INLINE ulong
_sc_submul(nn_ptr X, nn_srcptr P, slong N, ulong q)
{
    ulong cy = 0, hi, lo, b;
    slong i;

    for (i = 0; i < N; i++)
    {
        umul_ppmm(hi, lo, P[i], q);
        lo += cy;
        hi += (lo < cy);
        b = (X[i] < lo);
        X[i] -= lo;
        cy = hi + b;
    }
    return cy;
}

/* one integral limb (x->exp == 1), exact, n + 1 <= SC_PI2_LIMBS:
   q from (x1, x0) (R1, R0) without the low products, a lower
   approximation of x 2/pi within 5/B; P = pi/2 at N = n + 1 fraction
   limbs, truncated, t = X - q P in [0, pi/2 (1 + 5/B)).  Errors:
   (q + 1)(pi/2 - P) < B^-n, x truncated below B^-N, v truncated to n
   limbs: under 3 ulps of B^-n. */
static void
_sc_reduce1(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong n)
{
    ulong X[SC_PI2_LIMBS + 1], V[SC_PI2_LIMBS + 1];
    nn_srcptr P;
    slong N = n + 1;
    ulong q, h11, l11, h10, l10, h01, l01, c, u, top, b;
    int a, s = 1, fold, xneg = x->negative;

    _mp_real_elem_copy(X, N + 1, N, x);
    P = _sc_pi2_frac + SC_PI2_LIMBS - N;

    umul_ppmm(h11, l11, X[N], SC_2_DIV_PI_1);
    umul_ppmm(h10, l10, X[N], SC_2_DIV_PI_0);
    umul_ppmm(h01, l01, X[N - 1], SC_2_DIV_PI_1);
    (void) l10;
    (void) l01;
    add_ssaaaa(c, u, 0, l11, 0, h10);
    add_ssaaaa(c, u, c, u, 0, h01);
    q = h11 + c;
    a = (int) (q & 3);

    /* t = X - q (1 + P B^-N) */
    switch (N)
    {
        case 2: b = _sc_submul(X, P, 2, q); break;
        case 3: b = _sc_submul(X, P, 3, q); break;
        case 4: b = _sc_submul(X, P, 4, q); break;
        case 5: b = _sc_submul(X, P, 5, q); break;
        default: b = mpn_submul_1(X, P, N, q); break;
    }
    X[N] = X[N] - q - b;
    FLINT_ASSERT(X[N] <= 1);

    /* fold when t is above about pi/4: either way v <= pi/4 + 1/B */
    fold = (X[N] != 0) || (X[N - 1] > SC_PI4_TOP);

    if (fold)
    {
        /* v = (1 + P B^-N) - t: negative (and tiny) only when q was
           one low */
        b = mpn_sub_n(V, P, X, N);
        top = UWORD(1) - X[N] - b;
        V[N] = top;
        if ((slong) top < 0)
        {
            mpn_neg(V, V, N + 1);
            s = 1;
        }
        else
            s = -1;
        a = (a + 1) & 3;
        FLINT_ASSERT(V[N] == 0);
        _sc_eval(rs, rc, V + 1, n, 3, a, s, xneg);
    }
    else
        _sc_eval(rs, rc, X + 1, n, 3, a, s, xneg);
}

#endif

/* |x| >= 1 (m = exp >= 1 integral limbs), exact: reduce mod pi/2 and
   evaluate at n fraction limbs */
static void
_sc_reduce(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong n)
{
    slong m = x->exp, N = m + n + 1, xl;
    nn_ptr X, P, V, q;
    int a, s = 1, xneg = x->negative;
    TMP_INIT;

#if FLINT_BITS == 64
    if (m == 1 && n + 1 <= SC_PI2_LIMBS)
    {
        _sc_reduce1(rs, rc, x, n);
        return;
    }
#endif

    TMP_START;
    X = TMP_ALLOC((m + N + 1) * sizeof(ulong));
    P = TMP_ALLOC((N + 1) * sizeof(ulong));
    V = TMP_ALLOC((N + 1) * sizeof(ulong));
    q = TMP_ALLOC((2 * m + 4) * sizeof(ulong));

    /* X = |x| with m integral and N fraction limbs, truncated below
       B^-N (an error below one ulp of B^-N) */
    xl = m + N;
    _mp_real_elem_copy(X, xl, N, x);

    /* P = 2 floor(pi/4 B^N) <= pi/2 B^N, below by less than 2 */
    P[N] = mpn_lshift(P, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, N), N, 1);

    /* q = floor of a lower approximation of x 2/pi within 2/B:
       the top m + 1 limbs of X (m integral, one fraction) times
       floor(2/pi B^(m + 2)), limbs m + 3 up of the product */
    {
        slong k = m + 2;
        nn_ptr pr = TMP_ALLOC((2 * m + 3) * sizeof(ulong));

        flint_mpn_mul(pr, _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, k), k,
            X + N - 1, m + 1);
        flint_mpn_copyi(q, pr + m + 3, m);
    }

    /* a = q mod 4; t = X - q P, in [0, pi/2 (1 + 2/B)): the low N + 1
       limbs of X, the rest zero */
    a = (int) (q[0] & 3);
    {
        slong ql = m;
        nn_ptr T = TMP_ALLOC((ql + N + 1) * sizeof(ulong));

        while (ql > 0 && q[ql - 1] == 0)
            ql--;
        if (ql == 1)
        {
            ulong cy = mpn_submul_1(X, P, N + 1, q[0]);
            FLINT_ASSERT(cy == 0 || xl > N + 1);
            if (xl > N + 1)
                mpn_sub_1(X + N + 1, X + N + 1, xl - N - 1, cy);
        }
        else if (ql > 0)
        {
            if (ql >= N + 1)
                flint_mpn_mul(T, q, ql, P, N + 1);
            else
                flint_mpn_mul(T, P, N + 1, q, ql);
            /* q P <= X, so T has at most xl significant limbs */
            FLINT_ASSERT(ql + N + 1 <= xl || T[xl] == 0);
            mpn_sub(X, X, xl, T, FLINT_MIN(xl, ql + N + 1));
        }
    }
    FLINT_ASSERT(X[N] <= 1);

    /* fold: t > pi/4 (pi/4 < t when the units limb is set, else by the
       fraction limbs against those of pi/4 = P/2) */
    {
        int fold = (X[N] != 0);

        if (!fold)
        {
            nn_ptr H = V;
            mpn_rshift(H, P, N + 1, 1);
            fold = (mpn_cmp(X, H, N) > 0);
        }

        if (fold)
        {
            /* v = P - t: negative (and tiny) only when q was one low */
            if (mpn_sub_n(V, P, X, N + 1))
            {
                mpn_neg(V, V, N + 1);
                s = 1;
            }
            else
                s = -1;
            a = (a + 1) & 3;
        }
        else
            flint_mpn_copyi(V, X, N + 1);
    }

    /* v in [0, pi/4] at n fraction limbs: the reduction error
       q (pi/2 - P) (twice for the fold) < 2 B^(m - N) = 2 B^-(n+1), the
       truncations of x and v: 3 ulps of B^-n in all */
    _sc_eval(rs, rc, V + (N - n), n, 3, a, s, xneg);
    TMP_END;
}

void
mp_real_sin_cos_bits(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    slong e;
    ulong xerr = x->err;
    slong xanc = x->exp - x->size;

    prec = FLINT_MAX(prec, 2);

    /* x = 0 exactly */
    if (x->size == 0 && x->err == 0)
    {
        if (rs != NULL)
            mp_real_zero(rs);
        if (rc != NULL)
            mp_real_set_ui(rc, 1);
        return;
    }

    /* |x| < 2^e over the whole ball */
    e = mp_real_abs_bound_lt_2exp_si(x);

    /* a ball around zero, or one too wide for a relative accuracy of
       two bits within (-1, 1): sin in [+-min(2^e, 1)] (|sin x| <= |x|),
       cos = 1 +- min(2^(2e-1), 2) */
    if (x->size == 0 || (x->err != 0 && e <= 0
            && mp_real_rel_radius_lt_2exp_si(x) > -2))
    {
        if (rs != NULL)
        {
            mp_real_zero(rs);
            mp_real_add_error_2exp_si(rs, FLINT_MIN(e, 0));
        }
        if (rc != NULL)
        {
            mp_real_set_ui(rc, 1);
            mp_real_add_error_2exp_si(rc, FLINT_MIN(2 * e - 1, 1));
        }
        return;
    }

    /* ridiculously large arguments */
    if (FLINT_BITS * (x->exp - 1) >= FLINT_MAX(65536, 4 * prec))
    {
        if (rs != NULL)
        {
            mp_real_zero(rs);
            mp_real_add_error_2exp_si(rs, 0);
        }
        if (rc != NULL)
        {
            mp_real_zero(rc);
            mp_real_add_error_2exp_si(rc, 0);
        }
        return;
    }

    /* the midpoint, exactly */
    mid = *x;
    mid.err = 0;

    /* an inexact x: the absolute radius 2^-acc bounds the output error;
       for |x| < 1 the sine's relative accuracy is that of x */
    if (x->err != 0)
    {
        slong rel = mp_real_rel_radius_lt_2exp_si(x);   /* rad < 2^rel |x| */
        slong acc;

        if (x->exp <= 0)
            acc = -rel;                                 /* relative, for sin */
        else
            acc = -(rel + FLINT_BITS * x->exp);         /* absolute */

        if (acc < 2)
        {
            /* no useful information */
            if (rs != NULL)
            {
                mp_real_zero(rs);
                mp_real_add_error_2exp_si(rs, 0);
            }
            if (rc != NULL)
            {
                mp_real_zero(rc);
                mp_real_add_error_2exp_si(rc, 0);
            }
            return;
        }

        prec = FLINT_MIN(prec, acc + 8);
    }

    /* (the outputs may alias x: from here on, only mid, xerr and xanc
       are read, and each path copies mid before writing an output) */
    if (x->exp <= 0)
    {
        /* |x| < 2^emid, emid <= 0 */
        slong emid = FLINT_BITS * (x->exp - 1) + FLINT_BIT_COUNT(x->d[x->size - 1]);
        slong z = -emid;

        if (2 * z >= prec + 3)
            _sc_small1(rs, rc, &mid, emid);
        else if (4 * z >= prec + 3)
            _sc_small2(rs, rc, &mid, emid,
                (prec + z + 6 + FLINT_BITS - 1) / FLINT_BITS,
                mp_real_prec_bits(prec + 4));
        else
            _sc_unit(rs, rc, &mid, _sc_limbs(prec + z));
    }
    else
    {
        _sc_reduce(rs, rc, &mid, _sc_limbs(prec));
    }

    /* the radius of x */
    if (xerr != 0)
    {
        if (rs != NULL)
            _mp_real_elem_add_rad(rs, xerr, xanc);
        if (rc != NULL)
            _mp_real_elem_add_rad(rc, xerr, xanc);
    }
}
