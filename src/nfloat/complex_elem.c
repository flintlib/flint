/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "mpn_extras.h"
#include "nfloat.h"

/*
    Complex elementary functions, as compositions of the real functions
    using formulas without catastrophic cancellation except where the
    function itself is ill-conditioned (no internal guard bits apart from
    the exact evaluation of |z|^2 - 1 near the unit circle; each
    component typically has an error of a few ulp relative to the
    magnitude of the result, but small components of large results can
    have large relative errors, and directed rounding is not respected).
    Real arguments in the real domain of the function use the real version.

    Status flags of the real operations are or'ed together rather than
    checked after each operation: on over- or underflow (without
    NFLOAT_ALLOW_INF or NFLOAT_ALLOW_UNDERFLOW) the result is garbage
    and GR_UNABLE is returned. Squaring intermediate quantities reduces
    the usable exponent range by a bit, which does not matter in practice.
*/

#define RE(z) NFLOAT_COMPLEX_RE(z, ctx)
#define IM(z) NFLOAT_COMPLEX_IM(z, ctx)
#define CTMP(T) ulong T[2 * NFLOAT_MAX_ALLOC]
#define RTMP(T) ulong T[NFLOAT_MAX_ALLOC]

/* Temporaries are initialized to zero, so that values left unwritten by
   failing operations remain valid input for subsequent operations
   (whose results are then discarded). */
#define RZ(T) (NFLOAT_EXP(T) = NFLOAT_EXP_ZERO, NFLOAT_SGNBIT(T) = 0)
#define CZ(T) (RZ(RE(T)), RZ(IM(T)))

/* x = 2^k x in place */
static inline int
_nfloat_mul_2exp_ip(nfloat_ptr x, slong k, gr_ctx_t ctx)
{
    if (!NFLOAT_IS_SPECIAL(x))
    {
        NFLOAT_EXP(x) += k;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(x, ctx);
    }
    return GR_SUCCESS;
}

typedef enum
{
    C_EXP, C_EXPM1, C_EXP_PI_I, C_EXP2, C_EXP10,
    C_LOG, C_LOG1P, C_LOG_PI_I, C_LOG2, C_LOG10,
    C_SIN, C_COS, C_TAN, C_COT, C_SEC, C_CSC,
    C_SIN_PI, C_COS_PI, C_TAN_PI, C_COT_PI, C_SEC_PI, C_CSC_PI,
    C_SINC, C_SINC_PI,
    C_SINH, C_COSH, C_TANH, C_COTH, C_SECH, C_CSCH,
    C_ASIN, C_ACOS, C_ATAN, C_ACOT, C_ASEC, C_ACSC,
    C_ASINH, C_ACOSH, C_ATANH, C_ACOTH, C_ASECH, C_ACSCH,
    C_ASIN_PI, C_ACOS_PI, C_ATAN_PI, C_ACOT_PI, C_ASEC_PI, C_ACSC_PI
}
cfunc_t;

/* real domains */
enum { D_NONE, D_ALL, D_POS, D_GT_M1, D_ABS_LE1, D_ABS_LT1, D_ABS_GE1, D_ABS_GT1, D_GE1, D_ASECH };

typedef int (*_real_func_t)(nfloat_ptr, nfloat_srcptr, gr_ctx_t);

static const struct
{
    _real_func_t real;
    unsigned char domain;
}
_cfunc_real[] =
{
    { nfloat_exp, D_ALL }, { nfloat_expm1, D_ALL }, { NULL, D_NONE }, { nfloat_exp2, D_ALL }, { nfloat_exp10, D_ALL },
    { nfloat_log, D_POS }, { nfloat_log1p, D_GT_M1 }, { NULL, D_NONE }, { nfloat_log2, D_POS }, { nfloat_log10, D_POS },
    { nfloat_sin, D_ALL }, { nfloat_cos, D_ALL }, { nfloat_tan, D_ALL }, { nfloat_cot, D_ALL }, { nfloat_sec, D_ALL }, { nfloat_csc, D_ALL },
    { nfloat_sin_pi, D_ALL }, { nfloat_cos_pi, D_ALL }, { nfloat_tan_pi, D_ALL }, { nfloat_cot_pi, D_ALL }, { nfloat_sec_pi, D_ALL }, { nfloat_csc_pi, D_ALL },
    { nfloat_sinc, D_ALL }, { nfloat_sinc_pi, D_ALL },
    { nfloat_sinh, D_ALL }, { nfloat_cosh, D_ALL }, { nfloat_tanh, D_ALL }, { nfloat_coth, D_ALL }, { nfloat_sech, D_ALL }, { nfloat_csch, D_ALL },
    { nfloat_asin, D_ABS_LE1 }, { nfloat_acos, D_ABS_LE1 }, { nfloat_atan, D_ALL }, { nfloat_acot, D_ALL }, { nfloat_asec, D_ABS_GE1 }, { nfloat_acsc, D_ABS_GE1 },
    { nfloat_asinh, D_ALL }, { nfloat_acosh, D_GE1 }, { nfloat_atanh, D_ABS_LT1 }, { nfloat_acoth, D_ABS_GT1 }, { nfloat_asech, D_ASECH }, { nfloat_acsch, D_ALL },
    { nfloat_asin_pi, D_ABS_LE1 }, { nfloat_acos_pi, D_ABS_LE1 }, { nfloat_atan_pi, D_ALL }, { nfloat_acot_pi, D_ALL }, { nfloat_asec_pi, D_ABS_GE1 }, { nfloat_acsc_pi, D_ABS_GE1 },
};

/* whether the real number a lies in the real domain d */
static int
_in_domain(nfloat_srcptr a, int d, gr_ctx_t ctx)
{
    int c, s;
    RTMP(one);

    if (d == D_ALL)
        return 1;
    if (d == D_NONE)
        return 0;

    s = NFLOAT_IS_ZERO(a) ? 0 : (NFLOAT_SGNBIT(a) ? -1 : 1);

    if (d == D_POS)
        return s > 0;

    nfloat_one(one, ctx);
    if (d == D_GT_M1)
    {
        if (s >= 0)
            return 1;
        nfloat_cmpabs(&c, a, one, ctx);
        return c < 0;
    }

    nfloat_cmpabs(&c, a, one, ctx);

    switch (d)
    {
        case D_ABS_LE1: return c <= 0;
        case D_ABS_LT1: return c < 0;
        case D_ABS_GE1: return c >= 0;
        case D_ABS_GT1: return c > 0;
        case D_GE1: return s > 0 && c >= 0;
        default: return s > 0 && c <= 0;
    }
}

/* res = i^k z */
static void
_c_mul_i_pow(nfloat_complex_ptr res, nfloat_complex_srcptr z, int k, gr_ctx_t ctx)
{
    CTMP(t);

    k &= 3;
    if (k == 0)
    {
        nfloat_complex_set(res, (nfloat_complex_ptr) z, ctx);
        return;
    }
    if (k == 2)
    {
        nfloat_complex_neg(res, z, ctx);
        return;
    }

    /* i (a + bi) = -b + ai, -i (a + bi) = b - ai */
    nfloat_set(RE(t), IM(z), ctx);
    nfloat_set(IM(t), RE(z), ctx);
    if (k == 1)
        nfloat_neg(RE(t), RE(t), ctx);
    else
        nfloat_neg(IM(t), IM(t), ctx);
    nfloat_complex_set(res, t, ctx);
}

/* res = x c (c real) */
static int
_c_mul_real(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_srcptr c, gr_ctx_t ctx)
{
    int status;
    status = nfloat_mul(RE(res), RE(x), c, ctx);
    status |= nfloat_mul(IM(res), IM(x), c, ctx);
    return status;
}

static int
_c_div_real(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_srcptr c, gr_ctx_t ctx)
{
    int status;
    status = nfloat_div(RE(res), RE(x), c, ctx);
    status |= nfloat_div(IM(res), IM(x), c, ctx);
    return status;
}

/* sin and cos of a + bi (of pi (a + bi) if pi); s or c may be NULL */
static int
_c_sin_cos(nfloat_complex_ptr s, nfloat_complex_ptr c, nfloat_srcptr a, nfloat_srcptr b, int pi, gr_ctx_t ctx)
{
    RTMP(sa); RTMP(ca); RTMP(sb); RTMP(cb); RTMP(t);
    int status = GR_SUCCESS;

    RZ(sa); RZ(ca); RZ(sb); RZ(cb); RZ(t);

    /* sin(a + bi) = sin a cosh b + i cos a sinh b,
       cos(a + bi) = cos a cosh b - i sin a sinh b */
    if (pi)
    {
        status |= nfloat_sin_cos_pi(sa, ca, a, ctx);
        status |= nfloat_pi(t, ctx);
        status |= nfloat_mul(t, t, b, ctx);
        status |= nfloat_sinh_cosh(sb, cb, t, ctx);
    }
    else
    {
        status |= nfloat_sin_cos(sa, ca, a, ctx);
        status |= nfloat_sinh_cosh(sb, cb, b, ctx);
    }

    if (s != NULL)
    {
        status |= nfloat_mul(RE(s), sa, cb, ctx);
        status |= nfloat_mul(IM(s), ca, sb, ctx);
    }

    if (c != NULL)
    {
        status |= nfloat_mul(RE(c), ca, cb, ctx);
        status |= nfloat_mul(IM(c), sa, sb, ctx);
        nfloat_neg(IM(c), IM(c), ctx);
    }

    return status;
}

/* sec(a + bi) or csc(a + bi) (of pi (a + bi) if pi) */
static int
_c_sec_csc(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, int pi, int csc, gr_ctx_t ctx)
{
    RTMP(sa); RTMP(ca); RTMP(sb); RTMP(cb); RTMP(t);
    nfloat_ptr x, y;
    int status = GR_SUCCESS;

    RZ(sa); RZ(ca); RZ(sb); RZ(cb); RZ(t);

    /* sec(a + bi) = (cos a cosh b + i sin a sinh b) / (cos^2 a + sinh^2 b),
       csc(a + bi) = (sin a cosh b - i cos a sinh b) / (sin^2 a + sinh^2 b),
       with sums of nonnegative terms in the denominators */
    if (pi)
    {
        status |= nfloat_sin_cos_pi(sa, ca, a, ctx);
        status |= nfloat_pi(t, ctx);
        status |= nfloat_mul(t, t, b, ctx);
        status |= nfloat_sinh_cosh(sb, cb, t, ctx);
    }
    else
    {
        status |= nfloat_sin_cos(sa, ca, a, ctx);
        status |= nfloat_sinh_cosh(sb, cb, b, ctx);
    }

    /* (x, y) = (cos a, sin a) or (sin a, -cos a) */
    x = ca;
    y = sa;
    if (csc)
    {
        x = sa;
        y = ca;
        nfloat_neg(y, y, ctx);
    }

    status |= nfloat_sqr(t, x, ctx);
    status |= nfloat_mul(x, x, cb, ctx);
    status |= nfloat_sqr(cb, sb, ctx);
    status |= nfloat_add(t, t, cb, ctx);
    status |= nfloat_mul(y, y, sb, ctx);
    status |= nfloat_div(RE(res), x, t, ctx);
    status |= nfloat_div(IM(res), y, t, ctx);
    return status;
}

/* tan(a + bi) or cot(a + bi) (of pi (a + bi) if pi) */
static int
_c_tan_cot(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, int pi, int cot, gr_ctx_t ctx)
{
    RTMP(s); RTMP(c); RTMP(t); RTMP(u); RTMP(v);
    int status = GR_SUCCESS;

    RZ(s); RZ(c); RZ(t); RZ(u); RZ(v);

    if (pi)
    {
        status |= nfloat_pi(v, ctx);
        status |= nfloat_mul(v, v, b, ctx);
        b = v;
    }

    /* tan(bi) = i tanh(b), cot(bi) = -i coth(b) */
    if (NFLOAT_IS_ZERO(a))
    {
        nfloat_zero(RE(res), ctx);
        if (cot)
        {
            status |= nfloat_coth(IM(res), b, ctx);
            nfloat_neg(IM(res), IM(res), ctx);
        }
        else
        {
            status |= nfloat_tanh(IM(res), b, ctx);
        }
        return status;
    }

    if (!NFLOAT_IS_ZERO(b) && NFLOAT_EXP(b) >= 0)
    {
        int neg = NFLOAT_SGNBIT(b) ^ cot;

        /* |b| >= 1/2: with E = exp(-2|b|) <= exp(-1) and s, c the sine and
           cosine of 2a,
               tan(a + bi) = (2E s + i sgn(b) (1 - E^2)) / (1 + 2E c + E^2),
               cot(a + bi) = (2E s - i sgn(b) (1 - E^2)) / (1 - 2E c + E^2),
           where the denominators are >= (1 - E)^2 > 0.4. All quantities
           are bounded (the formula below would overflow in sinh^2 b long
           before the real part underflows), and the imaginary part is
           computed as 1 plus a correction of order E. */
        status |= nfloat_mul_2exp_si(t, a, 1, ctx);
        if (pi)
            status |= nfloat_sin_cos_pi(s, c, t, ctx);
        else
            status |= nfloat_sin_cos(s, c, t, ctx);

        status |= nfloat_mul_2exp_si(t, b, 1, ctx);
        NFLOAT_SGNBIT(t) = 1;
        if (nfloat_exp(t, t, ctx) != GR_SUCCESS)
        {
            /* E underflows: so does the real part, unless s = 0 */
            nfloat_zero(t, ctx);
            if (!NFLOAT_IS_ZERO(s))
                status |= GR_UNABLE;
        }
        /* (E^2 is negligible, and may underflow when E does not) */
        if (NFLOAT_EXP(t) > -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx) - 2)
            status |= nfloat_sqr(u, t, ctx);
        else
            nfloat_zero(u, ctx);
        status |= _nfloat_mul_2exp_ip(t, 1, ctx);
        status |= nfloat_mul(s, s, t, ctx);
        status |= nfloat_mul(c, c, t, ctx);
        if (cot)
            nfloat_neg(c, c, ctx);
        /* with c' = +-c: D = 1 + 2E c' + E^2, and
           (1 - E^2) / D = 1 - (2E c' + 2E^2) / D is computed in the
           second form (more accurate for small E) */
        nfloat_one(t, ctx);
        status |= nfloat_add(v, t, u, ctx);
        status |= nfloat_add(v, v, c, ctx);
        status |= nfloat_div(RE(res), s, v, ctx);
        status |= _nfloat_mul_2exp_ip(u, 1, ctx);
        status |= nfloat_add(u, u, c, ctx);
        status |= nfloat_div(u, u, v, ctx);
        status |= nfloat_sub(IM(res), t, u, ctx);
        NFLOAT_SGNBIT(IM(res)) = neg;
    }
    else
    {
        /* |b| < 1/2:
               tan(a + bi) = (sin a cos a + i sinh b cosh b) / (cos^2 a + sinh^2 b),
               cot(a + bi) = (sin a cos a - i sinh b cosh b) / (sin^2 a + sinh^2 b),
           with sums of nonnegative terms in the denominators */
        if (pi)
            status |= nfloat_sin_cos_pi(s, c, a, ctx);
        else
            status |= nfloat_sin_cos(s, c, a, ctx);
        status |= nfloat_sinh_cosh(t, u, b, ctx);
        status |= nfloat_mul(u, u, t, ctx);
        status |= nfloat_sqr(t, t, ctx);
        status |= nfloat_mul(v, s, c, ctx);
        status |= nfloat_sqr(c, cot ? s : c, ctx);
        status |= nfloat_add(c, c, t, ctx);
        status |= nfloat_div(RE(res), v, c, ctx);
        status |= nfloat_div(IM(res), u, c, ctx);
        if (cot)
            nfloat_neg(IM(res), IM(res), ctx);
    }

    return status;
}

/* exp(a + bi) */
static int
_c_exp(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    RTMP(E); RTMP(s); RTMP(c);
    int status;

    RZ(E); RZ(s); RZ(c);

    status = nfloat_sin_cos(s, c, b, ctx);
    if (NFLOAT_IS_ZERO(a))
    {
        nfloat_set(RE(res), c, ctx);
        nfloat_set(IM(res), s, ctx);
        return status;
    }
    status |= nfloat_exp(E, a, ctx);
    status |= nfloat_mul(RE(res), E, c, ctx);
    status |= nfloat_mul(IM(res), E, s, ctx);
    return status;
}

/* y = x with the precision increased from n to w limbs (exactly) */
static void
_nfloat_widen(nfloat_ptr y, nfloat_srcptr x, slong n, slong w)
{
    NFLOAT_EXP(y) = NFLOAT_EXP(x);
    NFLOAT_SGNBIT(y) = NFLOAT_SGNBIT(x);
    if (!NFLOAT_IS_SPECIAL(x))
    {
        flint_mpn_zero(NFLOAT_D(y), w - n);
        flint_mpn_copyi(NFLOAT_D(y) + w - n, NFLOAT_D(x), n);
    }
}

/* y = x with the precision reduced from w to n limbs (by truncation) */
static void
_nfloat_narrow(nfloat_ptr y, nfloat_srcptr x, slong n, slong w)
{
    NFLOAT_EXP(y) = NFLOAT_EXP(x);
    NFLOAT_SGNBIT(y) = NFLOAT_SGNBIT(x);
    if (!NFLOAT_IS_SPECIAL(x))
        flint_mpn_copyi(NFLOAT_D(y), NFLOAT_D(x) + w - n, n);
}

/*
    d = p^2 + q^2 - 1 where 1/2 <= |p| < 2 and q = 0 or EXP(q) <= EXP(p),
    computed exactly in fixed point (with 2^(-192n) as unit and an
    integral limb on top) and then truncated to n limbs. Returns 0 without
    output if q^2 lies below the fixed-point range (|q| < 2^(-32n)).
*/
static int
_c_abs2_m1_exact(nfloat_ptr d, nfloat_srcptr p, nfloat_srcptr q, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong W = 3 * n + 1, top, e, sh, off;
    ulong Xbuf[4 * NFLOAT_MAX_LIMBS + 1];
    ulong T[2 * NFLOAT_MAX_LIMBS + 1];
    ulong * X = Xbuf + n;
    unsigned int bs, c;
    int neg;

    if (!NFLOAT_IS_ZERO(q))
    {
        sh = FLINT_BITS * n + 2 * NFLOAT_EXP(q);
        if (sh < 0)
            return 0;
    }
    else
        sh = 0;

    flint_mpn_zero(Xbuf, 2 * n);

    /* p^2 - 1: p^2 < 4, so the top limb is in [0, 3] before subtracting */
    flint_mpn_sqr(T, NFLOAT_D(p), n);
    if (NFLOAT_EXP(p) == 1)
        X[3 * n] = mpn_lshift(X + n, T, 2 * n, 2);
    else
    {
        flint_mpn_copyi(X + n, T, 2 * n);
        X[3 * n] = 0;
    }
    X[3 * n] -= 1;

    /* + q^2 (two's complement) */
    if (!NFLOAT_IS_ZERO(q))
    {
        ulong cy;

        off = sh / FLINT_BITS;
        bs = sh % FLINT_BITS;
        flint_mpn_sqr(T, NFLOAT_D(q), n);
        T[2 * n] = 0;
        if (bs != 0)
            mpn_lshift(T, T, 2 * n + 1, bs);
        cy = mpn_add_n(X + off, X + off, T, 2 * n + 1);
        if (off + 2 * n + 1 < W)
            mpn_add_1(X + off + 2 * n + 1, X + off + 2 * n + 1, W - off - 2 * n - 1, cy);
    }

    neg = X[W - 1] >> (FLINT_BITS - 1);
    if (neg)
        mpn_neg(X, X, W);

    top = W - 1;
    while (top >= 0 && X[top] == 0)
        top--;

    if (top < 0)
    {
        nfloat_zero(d, ctx);
        return 1;
    }

    /* the top n limbs, from X[top - n], ..., X[top] (padded with zeros;
       addressed from Xbuf = X - n, which gcc can bound) */
    c = flint_clz(X[top]);
    e = FLINT_BITS * (top + 1) - c - 3 * FLINT_BITS * n;
    if (c != 0)
    {
        mpn_lshift(T, Xbuf + top, n + 1, c);
        flint_mpn_copyi(NFLOAT_D(d), T + 1, n);
    }
    else
        flint_mpn_copyi(NFLOAT_D(d), Xbuf + top + 1, n);

    NFLOAT_EXP(d) = e;
    NFLOAT_SGNBIT(d) = neg;
    return 1;
}

/*
    d = |z|^2 - 1, with an error small relative to |d|, where either
    (if log1p) z = 1 + p + qi with |p| < 1/2 and |q| < 1, or z = p + qi
    with 1/2 <= |p| < 2 and q = 0 or EXP(q) <= EXP(p).

    The second case is done exactly by the fixed-point code above unless
    |q| is tiny. Otherwise, the direct formula d = u + v with
    u = p (2 + p) respectively (|p| - 1)(|p| + 1) (with |p| - 1 exact)
    and v = q^2 has a small relative error unless u < 0 and there is
    cancellation. In that case the terms that can cancel are instead added
    exactly at 2n + 1 limbs: p^2 - 1 is exact, and in log1p mode 2p + q^2
    is exact unless the exponents differ by more than a limb (when there
    is no cancellation). The remaining term is then added with a single
    rounding.
*/
static int
_c_abs2_m1(nfloat_ptr d, nfloat_srcptr p, nfloat_srcptr q, int log1p, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong w = 2 * n + 1;
    RTMP(U); RTMP(V);
    int status = GR_SUCCESS;

    if (!log1p && _c_abs2_m1_exact(d, p, q, ctx))
        return GR_SUCCESS;

    RZ(U); RZ(V);

    nfloat_one(V, ctx);
    if (log1p)
    {
        NFLOAT_EXP(V) = 2;
        status |= nfloat_add(U, V, p, ctx);
        status |= nfloat_mul(U, U, p, ctx);
    }
    else
    {
        RTMP(A);

        nfloat_abs(A, p, ctx);
        status |= nfloat_add(U, A, V, ctx);
        status |= nfloat_sub(A, A, V, ctx);
        status |= nfloat_mul(U, U, A, ctx);
    }
    status |= nfloat_sqr(V, q, ctx);
    status |= nfloat_add(d, U, V, ctx);

    /* no cancellation, or at most about one bit */
    if (!NFLOAT_SGNBIT(U) || NFLOAT_IS_ZERO(V) || w > NFLOAT_MAX_LIMBS ||
        (!NFLOAT_IS_ZERO(d) && NFLOAT_EXP(d) >= FLINT_MAX(NFLOAT_EXP(U), NFLOAT_EXP(V)) - 1))
        return status;

    {
        gr_ctx_t wctx;
        RTMP(P); RTMP(Q); RTMP(S); RTMP(T);

        RZ(S); RZ(T);

        nfloat_ctx_init(wctx, w * FLINT_BITS, NFLOAT_CTX_FLAGS(ctx));
        _nfloat_widen(P, p, n, w);
        _nfloat_widen(Q, q, n, w);

        if (log1p)
        {
            /* (2p + q^2) + p^2 */
            status |= nfloat_sqr(S, Q, wctx);
            status |= nfloat_mul_2exp_si(T, P, 1, wctx);
            status |= nfloat_add(S, S, T, wctx);
            status |= nfloat_sqr(T, P, wctx);
            status |= nfloat_add(S, S, T, wctx);
        }
        else
        {
            /* (p^2 - 1) + q^2 */
            status |= nfloat_sqr(S, P, wctx);
            nfloat_one(T, wctx);
            status |= nfloat_sub(S, S, T, wctx);
            status |= nfloat_sqr(T, Q, wctx);
            status |= nfloat_add(S, S, T, wctx);
        }

        _nfloat_narrow(d, S, n, w);
    }

    return status;
}

/* log(a + bi), or log(1 + a + bi) if log1p */
static int
_c_log(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, int log1p, gr_ctx_t ctx)
{
    RTMP(t); RTMP(u); RTMP(x1);
    nfloat_srcptr x, p, q;
    int status = GR_SUCCESS;

    RZ(t); RZ(u); RZ(x1);

    if (log1p)
    {
        nfloat_one(t, ctx);
        status |= nfloat_add(x1, t, a, ctx);
        x = x1;
    }
    else
    {
        x = a;
    }

    if (NFLOAT_IS_ZERO(x) && NFLOAT_IS_ZERO(b))
        return GR_DOMAIN;

    status |= nfloat_atan2(IM(res), b, x, ctx);

    /* Re = log(|z|^2) / 2 = log1p(|z|^2 - 1) / 2 */
    if (log1p && NFLOAT_EXP(a) <= -1 && NFLOAT_EXP(b) <= 0)
    {
        /* |a| < 1/2 (1 + a inexact), |b| < 1: |z|^2 - 1 in (-3/4, 9/4) */
        status |= _c_abs2_m1(t, a, b, 1, ctx);
        status |= nfloat_log1p(RE(res), t, ctx);
    }
    else
    {
        /* p = the larger of |x|, |b| (if log1p, 1 + a is exact unless
           |a| < 1/2 and |b| >= 1, when |z|^2 >= 5/4 is not near 1) */
        p = x; q = b;
        if (NFLOAT_IS_ZERO(x) || (!NFLOAT_IS_ZERO(b) && NFLOAT_EXP(b) > NFLOAT_EXP(x)))
        {
            p = b; q = x;
        }

        if (NFLOAT_EXP(p) == 0 || NFLOAT_EXP(p) == 1)
        {
            /* 1/2 <= |p| < 2, |z|^2 in [1/4, 8) */
            status |= _c_abs2_m1(t, p, q, 0, ctx);
            status |= nfloat_log1p(RE(res), t, ctx);
        }
        else
        {
            /* |z|^2 < 1/2 or |z|^2 >= 4 */
            status |= nfloat_sqr(t, x, ctx);
            status |= nfloat_sqr(u, b, ctx);
            status |= nfloat_add(t, t, u, ctx);
            status |= nfloat_log(RE(res), t, ctx);
        }
    }

    status |= _nfloat_mul_2exp_ip(RE(res), -1, ctx);
    return status;
}

/* atan(a + bi), b != 0, a + bi != +-i */
static int
_c_atan(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    RTMP(t); RTMP(u); RTMP(v); RTMP(one);
    nfloat_srcptr p, q;
    int status = GR_SUCCESS;

    RZ(t); RZ(u); RZ(v);

    /* atan is odd: b > 0 avoids cancellation in 1 + b below */
    if (NFLOAT_SGNBIT(b))
    {
        RTMP(na); RTMP(nb);

        nfloat_neg(na, a, ctx);
        nfloat_neg(nb, b, ctx);
        status = _c_atan(res, na, nb, ctx);
        nfloat_complex_neg(res, res, ctx);
        return status;
    }

    /* Re = atan2(2a, 1 - |z|^2) / 2,
       Im = log1p(4b / (a^2 + (1 - b)^2)) / 4, the argument of log1p
       positive */
    nfloat_one(one, ctx);
    status |= nfloat_sub(t, one, b, ctx);
    status |= nfloat_sqr(t, t, ctx);
    if (!NFLOAT_IS_ZERO(a))
    {
        status |= nfloat_sqr(v, a, ctx);
        status |= nfloat_add(t, t, v, ctx);
    }

    /* u = 1 - |z|^2, accurately (1 - b and 1 + b are only used for the
       cases without cancellation) */
    p = a; q = b;
    if (NFLOAT_IS_ZERO(a) || NFLOAT_EXP(b) > NFLOAT_EXP(a))
    {
        p = b; q = a;
    }
    if (NFLOAT_EXP(p) == 0 || NFLOAT_EXP(p) == 1)
    {
        status |= _c_abs2_m1(u, p, q, 0, ctx);
        nfloat_neg(u, u, ctx);
    }
    else
    {
        status |= nfloat_sub(u, one, b, ctx);
        status |= nfloat_add(v, one, b, ctx);
        status |= nfloat_mul(u, u, v, ctx);
        if (!NFLOAT_IS_ZERO(a))
        {
            status |= nfloat_sqr(v, a, ctx);
            status |= nfloat_sub(u, u, v, ctx);
        }
    }

    status |= nfloat_div(t, b, t, ctx);
    status |= _nfloat_mul_2exp_ip(t, 2, ctx);
    status |= nfloat_log1p(IM(res), t, ctx);
    status |= _nfloat_mul_2exp_ip(IM(res), -2, ctx);

    if (NFLOAT_IS_ZERO(a))
    {
        /* atan(bi) for |b| > 1 has real part sgn(b) pi/2 */
        if (NFLOAT_SGNBIT(u))
        {
            status |= nfloat_pi(RE(res), ctx);
            NFLOAT_EXP(RE(res)) -= 1;
        }
        else
            nfloat_zero(RE(res), ctx);
        return status;
    }

    status |= nfloat_mul_2exp_si(v, a, 1, ctx);
    status |= nfloat_atan2(RE(res), v, u, ctx);
    status |= _nfloat_mul_2exp_ip(RE(res), -1, ctx);
    return status;
}

/* square roots of 1 - z and 1 + z (or z - 1 and z + 1 if minus) */
static int
_c_sqrt_1pm(nfloat_complex_ptr s1, nfloat_complex_ptr s2, nfloat_complex_srcptr z, int minus, gr_ctx_t ctx)
{
    CTMP(t);
    int status;

    nfloat_complex_one(t, ctx);
    if (minus)
        status = nfloat_complex_sub(s1, z, t, ctx);
    else
        status = nfloat_complex_sub(s1, t, z, ctx);
    status |= nfloat_complex_add(s2, t, z, ctx);
    status |= nfloat_complex_sqrt(s1, s1, ctx);
    status |= nfloat_complex_sqrt(s2, s2, ctx);
    return status;
}

/* Kahan's formulas for asin, acos and acosh */
static int
_c_asin_acos(nfloat_complex_ptr res, nfloat_complex_srcptr z, cfunc_t f, gr_ctx_t ctx)
{
    CTMP(s1); CTMP(s2);
    RTMP(t); RTMP(u);
    int status;

    CZ(s1); CZ(s2);
    RZ(t); RZ(u);

    if (f == C_ACOSH)
    {
        /* s1 = sqrt(z - 1), s2 = sqrt(z + 1):
           Re = asinh(Re(conj(s1) s2)), Im = 2 atan2(Im s1, Re s2) */
        status = _c_sqrt_1pm(s1, s2, z, 1, ctx);
        status |= nfloat_mul(t, RE(s1), RE(s2), ctx);
        status |= nfloat_addmul(t, IM(s1), IM(s2), ctx);
        status |= nfloat_atan2(IM(res), IM(s1), RE(s2), ctx);
        status |= _nfloat_mul_2exp_ip(IM(res), 1, ctx);
        status |= nfloat_asinh(RE(res), t, ctx);
        return status;
    }

    /* s1 = sqrt(1 - z), s2 = sqrt(1 + z) */
    status = _c_sqrt_1pm(s1, s2, z, 0, ctx);

    if (f == C_ASIN)
    {
        /* Re = atan2(a, Re(s1 s2)), Im = asinh(Im(conj(s1) s2)) */
        status |= nfloat_mul(t, RE(s1), RE(s2), ctx);
        status |= nfloat_submul(t, IM(s1), IM(s2), ctx);
        status |= nfloat_mul(u, RE(s1), IM(s2), ctx);
        status |= nfloat_submul(u, IM(s1), RE(s2), ctx);
        status |= nfloat_atan2(RE(res), RE(z), t, ctx);
        status |= nfloat_asinh(IM(res), u, ctx);
    }
    else
    {
        /* Re = 2 atan2(Re s1, Re s2), Im = asinh(Im(conj(s2) s1)) */
        status |= nfloat_mul(u, RE(s2), IM(s1), ctx);
        status |= nfloat_submul(u, IM(s2), RE(s1), ctx);
        status |= nfloat_atan2(RE(res), RE(s1), RE(s2), ctx);
        status |= _nfloat_mul_2exp_ip(RE(res), 1, ctx);
        status |= nfloat_asinh(IM(res), u, ctx);
    }

    return status;
}

static int
_nfloat_complex_func(nfloat_complex_ptr res, nfloat_complex_srcptr x, cfunc_t f, gr_ctx_t ctx)
{
    CTMP(zbuf); CTMP(wbuf);
    RTMP(t);
    nfloat_complex_srcptr z;
    nfloat_complex_ptr w;
    nfloat_srcptr a = RE(x), b = IM(x);
    int recip_arg = 0, rot = 0, div_pi = 0;
    int status = GR_SUCCESS;

    CZ(wbuf);
    RZ(t);

    if ((NFLOAT_IS_SPECIAL(a) && !NFLOAT_IS_ZERO(a)) || (NFLOAT_IS_SPECIAL(b) && !NFLOAT_IS_ZERO(b)))
        return GR_UNABLE;

    /* real arguments in the real domain */
    if (NFLOAT_IS_ZERO(b) && _in_domain(a, _cfunc_real[f].domain, ctx))
    {
        status = _cfunc_real[f].real(RE(res), a, ctx);
        nfloat_zero(IM(res), ctx);
        return status;
    }

    /* reduce to fewer cases: functions of 1/z, hyperbolic functions as
       trigonometric functions of iz */
    if (f >= C_ASIN_PI)
    {
        div_pi = 1;
        f = (cfunc_t) (f - C_ASIN_PI + C_ASIN);
    }

    switch (f)
    {
        case C_ACOT: f = C_ATAN; recip_arg = 1; break;
        case C_ASEC: f = C_ACOS; recip_arg = 1; break;
        case C_ACSC: f = C_ASIN; recip_arg = 1; break;
        case C_ACOTH: f = C_ATANH; recip_arg = 1; break;
        case C_ASECH: f = C_ACOSH; recip_arg = 1; break;
        case C_ACSCH: f = C_ASINH; recip_arg = 1; break;
        default: break;
    }

    /* sinh z = -i sin(iz), cosh z = cos(iz), tanh z = -i tan(iz),
       coth z = i cot(iz), sech z = sec(iz), csch z = i csc(iz),
       asinh z = -i asin(iz), atanh z = -i atan(iz) */
    switch (f)
    {
        case C_SINH: f = C_SIN; rot = 3; break;
        case C_COSH: f = C_COS; rot = 4; break;
        case C_TANH: f = C_TAN; rot = 3; break;
        case C_COTH: f = C_COT; rot = 1; break;
        case C_SECH: f = C_SEC; rot = 4; break;
        case C_CSCH: f = C_CSC; rot = 1; break;
        case C_ASINH: f = C_ASIN; rot = 3; break;
        case C_ATANH: f = C_ATAN; rot = 3; break;
        default: break;
    }

    z = x;
    if (recip_arg)
    {
        status = nfloat_complex_inv(zbuf, x, ctx);
        if (status != GR_SUCCESS)
            return status;
        z = zbuf;
    }
    if (rot)
    {
        _c_mul_i_pow(zbuf, z, 1, ctx);
        z = zbuf;
    }

    /* write directly to res unless it aliases the input */
    w = (res == x) ? wbuf : res;

    a = RE(z);
    b = IM(z);

    if (NFLOAT_IS_ZERO(b) && f != C_EXP_PI_I && f != C_LOG_PI_I &&
        _in_domain(a, _cfunc_real[f].domain, ctx))
    {
        status |= _cfunc_real[f].real(RE(w), a, ctx);
        nfloat_zero(IM(w), ctx);
    }
    else switch (f)
    {
        case C_EXP:
            status |= _c_exp(w, a, b, ctx);
            break;

        case C_EXPM1:
            if (NFLOAT_EXP(a) <= -1 && NFLOAT_EXP(b) <= 0)
            {
                /* |a| < 1/2, |b| < 1:
                   Re = expm1(a) cos b - sin^2 b / (1 + cos b),
                   Im = (1 + expm1(a)) sin b */
                RTMP(s); RTMP(c); RTMP(u);

                RZ(s); RZ(c); RZ(u);

                status |= nfloat_sin_cos(s, c, b, ctx);
                status |= nfloat_expm1(t, a, ctx);
                nfloat_one(u, ctx);
                status |= nfloat_add(u, u, c, ctx);
                status |= nfloat_mul(RE(w), t, c, ctx);
                status |= nfloat_mul(c, s, s, ctx);
                status |= nfloat_div(c, c, u, ctx);
                status |= nfloat_sub(RE(w), RE(w), c, ctx);
                nfloat_one(u, ctx);
                status |= nfloat_add(t, t, u, ctx);
                status |= nfloat_mul(IM(w), s, t, ctx);
            }
            else
            {
                status |= _c_exp(w, a, b, ctx);
                nfloat_one(t, ctx);
                status |= nfloat_sub(RE(w), RE(w), t, ctx);
            }
            break;

        case C_EXP_PI_I:
            /* exp(pi i z) = exp(-pi b) (cos(pi a) + i sin(pi a)) */
            status |= nfloat_sin_cos_pi(IM(w), RE(w), a, ctx);
            if (!NFLOAT_IS_ZERO(b))
            {
                status |= nfloat_pi(t, ctx);
                status |= nfloat_mul(t, t, b, ctx);
                nfloat_neg(t, t, ctx);
                status |= nfloat_exp(t, t, ctx);
                status |= _c_mul_real(w, w, t, ctx);
            }
            break;

        case C_EXP2:
        case C_EXP10:
            status |= nfloat_set_ui(t, (f == C_EXP2) ? 2 : 10, ctx);
            status |= nfloat_log(t, t, ctx);
            status |= _c_mul_real(w, z, t, ctx);
            status |= _c_exp(w, RE(w), IM(w), ctx);
            break;

        case C_LOG:
        case C_LOG2:
        case C_LOG10:
        case C_LOG_PI_I:
            status |= _c_log(w, a, b, 0, ctx);
            if (f == C_LOG2 || f == C_LOG10)
            {
                status |= nfloat_set_ui(t, (f == C_LOG2) ? 2 : 10, ctx);
                status |= nfloat_log(t, t, ctx);
                status |= _c_div_real(w, w, t, ctx);
            }
            else if (f == C_LOG_PI_I)
            {
                /* log(z) / (pi i) */
                status |= nfloat_pi(t, ctx);
                status |= _c_div_real(w, w, t, ctx);
                _c_mul_i_pow(w, w, 3, ctx);
            }
            break;

        case C_LOG1P:
            status |= _c_log(w, a, b, 1, ctx);
            break;

        case C_SIN:
        case C_SIN_PI:
            status |= _c_sin_cos(w, NULL, a, b, f == C_SIN_PI, ctx);
            break;

        case C_COS:
        case C_COS_PI:
            status |= _c_sin_cos(NULL, w, a, b, f == C_COS_PI, ctx);
            break;

        case C_TAN:
        case C_TAN_PI:
            status |= _c_tan_cot(w, a, b, f == C_TAN_PI, 0, ctx);
            break;

        case C_COT:
        case C_COT_PI:
            status |= _c_tan_cot(w, a, b, f == C_COT_PI, 1, ctx);
            break;

        case C_SEC:
        case C_SEC_PI:
            status |= _c_sec_csc(w, a, b, f == C_SEC_PI, 0, ctx);
            break;

        case C_CSC:
        case C_CSC_PI:
            status |= _c_sec_csc(w, a, b, f == C_CSC_PI, 1, ctx);
            break;

        case C_SINC:
        case C_SINC_PI:
            status |= _c_sin_cos(w, NULL, a, b, f == C_SINC_PI, ctx);
            status |= nfloat_complex_div(w, w, z, ctx);
            if (f == C_SINC_PI)
            {
                status |= nfloat_pi(t, ctx);
                status |= _c_div_real(w, w, t, ctx);
            }
            break;

        case C_ASIN:
        case C_ACOS:
        case C_ACOSH:
            status |= _c_asin_acos(w, z, f, ctx);
            break;

        case C_ATAN:
            status |= _c_atan(w, a, b, ctx);
            break;

        default:
            return GR_UNABLE;
    }

    /* GR_DOMAIN from log at zero */
    if (status == GR_DOMAIN)
        return status;

    if (rot)
        _c_mul_i_pow(w, w, rot, ctx);

    if (div_pi)
    {
        status |= nfloat_pi(t, ctx);
        status |= _c_div_real(w, w, t, ctx);
    }

    if (w != res)
        nfloat_complex_set(res, w, ctx);

    /* (GR_DOMAIN from a division by a garbage zero is not meaningful) */
    return (status & GR_UNABLE) ? GR_UNABLE : status;
}

#define DEF(name, F) \
    int nfloat_complex_ ## name(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx) \
    { return _nfloat_complex_func(res, x, F, ctx); }

DEF(exp, C_EXP)
DEF(expm1, C_EXPM1)
DEF(exp_pi_i, C_EXP_PI_I)
DEF(exp2, C_EXP2)
DEF(exp10, C_EXP10)
DEF(log, C_LOG)
DEF(log1p, C_LOG1P)
DEF(log_pi_i, C_LOG_PI_I)
DEF(log2, C_LOG2)
DEF(log10, C_LOG10)
DEF(sin, C_SIN)
DEF(cos, C_COS)
DEF(tan, C_TAN)
DEF(cot, C_COT)
DEF(sec, C_SEC)
DEF(csc, C_CSC)
DEF(sin_pi, C_SIN_PI)
DEF(cos_pi, C_COS_PI)
DEF(tan_pi, C_TAN_PI)
DEF(cot_pi, C_COT_PI)
DEF(sec_pi, C_SEC_PI)
DEF(csc_pi, C_CSC_PI)
DEF(sinc, C_SINC)
DEF(sinc_pi, C_SINC_PI)
DEF(sinh, C_SINH)
DEF(cosh, C_COSH)
DEF(tanh, C_TANH)
DEF(coth, C_COTH)
DEF(sech, C_SECH)
DEF(csch, C_CSCH)
DEF(asin, C_ASIN)
DEF(acos, C_ACOS)
DEF(atan, C_ATAN)
DEF(acot, C_ACOT)
DEF(asec, C_ASEC)
DEF(acsc, C_ACSC)
DEF(asinh, C_ASINH)
DEF(acosh, C_ACOSH)
DEF(atanh, C_ATANH)
DEF(acoth, C_ACOTH)
DEF(asech, C_ASECH)
DEF(acsch, C_ACSCH)
DEF(asin_pi, C_ASIN_PI)
DEF(acos_pi, C_ACOS_PI)
DEF(atan_pi, C_ATAN_PI)
DEF(acot_pi, C_ACOT_PI)
DEF(asec_pi, C_ASEC_PI)
DEF(acsc_pi, C_ACSC_PI)

static int
_nfloat_complex_sin_cos(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, int pi, gr_ctx_t ctx)
{
    CTMP(s); CTMP(c);
    int status;

    CZ(s); CZ(c);

    if ((NFLOAT_IS_SPECIAL(RE(x)) && !NFLOAT_IS_ZERO(RE(x))) || (NFLOAT_IS_SPECIAL(IM(x)) && !NFLOAT_IS_ZERO(IM(x))))
        return GR_UNABLE;

    status = _c_sin_cos(s, c, RE(x), IM(x), pi, ctx);
    nfloat_complex_set(res1, s, ctx);
    nfloat_complex_set(res2, c, ctx);
    return status;
}

int
nfloat_complex_sin_cos(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_complex_sin_cos(res1, res2, x, 0, ctx);
}

int
nfloat_complex_sin_cos_pi(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_complex_sin_cos(res1, res2, x, 1, ctx);
}

int
nfloat_complex_sinh_cosh(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    CTMP(z);
    int status;

    /* sinh z = -i sin(iz), cosh z = cos(iz) */
    _c_mul_i_pow(z, x, 1, ctx);
    status = nfloat_complex_sin_cos(res1, res2, z, ctx);
    _c_mul_i_pow(res1, res1, 3, ctx);
    return status;
}

int
nfloat_complex_pow(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    CTMP(t);
    int status;

    CZ(t);

    if (NFLOAT_COMPLEX_IS_SPECIAL(x, ctx) || NFLOAT_COMPLEX_IS_SPECIAL(y, ctx))
    {
        if (NFLOAT_COMPLEX_IS_ZERO(y, ctx))
            return nfloat_complex_one(res, ctx);

        if (NFLOAT_COMPLEX_IS_ZERO(x, ctx))
        {
            /* 0^y = 0 for Re(y) > 0 */
            if (!NFLOAT_IS_SPECIAL(RE(y)) && !NFLOAT_SGNBIT(RE(y)))
                return nfloat_complex_zero(res, ctx);
            return GR_DOMAIN;
        }
    }

    /* real powers of positive reals */
    if (NFLOAT_IS_ZERO(IM(x)) && NFLOAT_IS_ZERO(IM(y)) && !NFLOAT_SGNBIT(RE(x)))
    {
        status = nfloat_pow(RE(res), RE(x), RE(y), ctx);
        nfloat_zero(IM(res), ctx);
        return status;
    }

    /* integer powers (|y| < 2^16) by repeated squaring */
    if (NFLOAT_IS_ZERO(IM(y)) && !NFLOAT_IS_SPECIAL(RE(y)) && NFLOAT_EXP(RE(y)) <= 16)
    {
        RTMP(u);

        if (nfloat_nint(u, RE(y), ctx) == GR_SUCCESS && nfloat_equal(u, RE(y), ctx) == T_TRUE)
        {
            fmpz_t e;

            fmpz_init(e);
            status = nfloat_get_fmpz(e, u, ctx);
            if (status == GR_SUCCESS)
                status = gr_generic_pow_fmpz(res, x, e, ctx);
            fmpz_clear(e);
            return status;
        }
    }

    /* exp(y log(x)) */
    status = nfloat_complex_log(t, x, ctx);
    status |= nfloat_complex_mul(t, t, y, ctx);
    status |= nfloat_complex_exp(res, t, ctx);
    return status;
}
