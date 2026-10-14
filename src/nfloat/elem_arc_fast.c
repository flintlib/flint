/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

#if FLINT_BITS == 64

/*
    In-register fast paths for the inverse trigonometric and hyperbolic
    functions at one and two limbs (64-bit only), at w = n + 1 fraction
    limbs, for the compositions of elem_misc.c without their inner
    context:

      asin, acos         atan2 of y = |x| and S = sqrt(1 - y^2)
      asec, acsc         atan2 of S = sqrt(y^2 - 1) and 1
      acot, atan_pi      atan2 of 1 and y
      asinh              log(y + sqrt(y^2 + 1)) for y >= 1/2, and
                         log1p(y + sqrt(1 + y^2) - 1) below
      acosh              log(y + sqrt(y^2 - 1))
      atanh              (log1p(y) - log(1 - y)) / 2
      acoth              (log(y + 1) - log(y - 1)) / 2 for y <= 2, else
                         atanh(1/y)
      asech              log1p(S) - log(y), S = sqrt(1 - y^2)
      acsch              log(1 + sqrt(1 + y^2)) - log(y) for y < 1, else
                         asinh(1/y)

    with the functions of pi taking c/4 + sigma atan(t) (2/pi) / 2 for
    c pi/4 + sigma atan(t). The two terms of a sum or difference of
    logarithms have the same sign, so nothing cancels, and the arguments
    of the logarithms near 1 are computed with absolute accuracy, which is
    relative accuracy for the result where it is that small.

    Quantities are kept as _nfloat_fxr_t: a normalized w-limb mantissa
    D in [1/2, 1), an exponent, and an error bound in units of B^-w of D.
    y^2 +- 1 is computed exactly; square roots by _nfloat_small_sqrt_2 / 3
    (which halve relative errors); quotients by one 4/2 resp. 6/3
    division. The kernels are those of mp_real (atan of a ratio, log1p)
    and elem.h (log), with error bounds ERR. Propagated errors are bounded
    by the derivative factors noted below, so that the final bound is a
    valid count of ulps of the fixed-point result.

    The paths decline (NFLOAT_FAST_DECLINE) when an intermediate would
    lose relative accuracy at this precision: results below 2^-32 resp.
    2^-40 where they are computed as small differences, |x| >= 2^62,
    and |x| < 2^-31 resp. 2^-39 (the tiny cases are mostly handled by the
    callers).
*/

typedef struct
{
    ulong d[3];
    slong e;
    ulong err;
}
_nfloat_fxr_struct;

typedef _nfloat_fxr_struct _nfloat_fxr_t[1];

#define FXR_MIN_E(n) (((n) == 1) ? -31 : -39)
#define FXR_MAX_SHIFT(n) (((n) == 1) ? 32 : 40)
#define FXR_MAX_E 62

#define ATAN_ERR(w) (((w) == 2) ? MP_REAL_SMALL_ATAN_ERR_2 : MP_REAL_SMALL_ATAN_ERR_3)
#define LOG1P_ERR(w) (((w) == 2) ? MP_REAL_SMALL_LOG_ERR_2 : MP_REAL_SMALL_LOG_ERR_3)
#define LOG_ERR(w) (((w) == 2) ? NFLOAT_SMALL_LOG_ERR_2 : NFLOAT_SMALL_LOG_ERR_3)
#define SQRT_ERR(w) (((w) == 2) ? NFLOAT_SMALL_SQRT_ERR_2 : NFLOAT_SMALL_SQRT_ERR_3)

/* r = |x| exactly */
FLINT_FORCE_INLINE void
_fxr_set_nfloat(_nfloat_fxr_t r, nfloat_srcptr x, slong n)
{
    r->d[0] = 0;
    r->d[1] = NFLOAT_D(x)[0];
    /* (d[2] is unused for n = 1; written anyway, as gcc cannot tell) */
    r->d[2] = (n == 2) ? NFLOAT_D(x)[1] : 0;
    r->e = NFLOAT_EXP(x);
    r->err = 0;
}

FLINT_FORCE_INLINE void
_fxr_one(_nfloat_fxr_t r, slong w)
{
    r->d[0] = 0;
    r->d[1] = 0;
    r->d[w - 1] = UWORD(1) << (FLINT_BITS - 1);
    r->e = 1;
    r->err = 0;
}

/* r = the integer (z, zn) times 2^e0, z nonzero, truncated to w limbs */
static void
_fxr_set_mpn(_nfloat_fxr_t r, nn_srcptr z, slong zn, slong e0, slong w)
{
    ulong t[8];
    unsigned int c;
    slong i;

    while (z[zn - 1] == 0)
        zn--;

    c = flint_clz(z[zn - 1]);
    for (i = 0; i <= w; i++)
        t[i] = (zn - 1 - w + i >= 0) ? z[zn - 1 - w + i] : 0;
    if (c != 0)
        mpn_lshift(t, t, w + 1, c);
    for (i = 0; i < w; i++)
        r->d[i] = t[i + 1];

    r->e = FLINT_BITS * zn - c + e0;
    r->err = 1;
}

/* r = 1 - y^2 (mode 0), y^2 - 1 (mode -1) or y^2 + 1 (mode 1) for y exact
   with n limbs, 0 <= K = FLINT_BITS 2n - 2e < FLINT_BITS (2n + 2) - 1;
   within one unit */
static void
_fxr_sqr_pm1(_nfloat_fxr_t r, const _nfloat_fxr_t y, int mode, slong n, slong w)
{
    ulong z[6], p[6];
    slong K = 2 * FLINT_BITS * n - 2 * y->e, zn = 2 * n + 2;

    flint_mpn_zero(z, zn);
    flint_mpn_zero(p, zn);
    if (n == 1)
        umul_ppmm(z[1], z[0], y->d[1], y->d[1]);
    else
        FLINT_MPN_SQR_2X2(z[3], z[2], z[1], z[0], y->d[2], y->d[1]);
    p[K / FLINT_BITS] = UWORD(1) << (K % FLINT_BITS);

    if (mode == 0)
        mpn_sub_n(z, p, z, zn);
    else if (mode < 0)
        mpn_sub_n(z, z, p, zn);
    else
        mpn_add_n(z, z, p, zn);

    _fxr_set_mpn(r, z, zn, 2 * y->e - 2 * FLINT_BITS * n, w);
    (void) w;
}

/* s = sqrt(u) */
FLINT_FORCE_INLINE void
_fxr_sqrt(_nfloat_fxr_t s, const _nfloat_fxr_t u, slong w)
{
    ulong t = u->e & 1;

    if (w == 2)
        _nfloat_small_sqrt_2(&s->d[1], &s->d[0], u->d[1], u->d[0], t);
    else
        _nfloat_small_sqrt_3(&s->d[2], &s->d[1], &s->d[0], u->d[2], u->d[1], u->d[0], t);

    /* the relative error of u halved: within u->err units of the
       mantissa in [1/2, 1) */
    s->e = (u->e + (slong) t) / 2;
    s->err = SQRT_ERR(w) + u->err;
}

/* (v, w) = floor(a.d B^w 2^-sh) */
FLINT_FORCE_INLINE void
_fxr_rshift_w(nn_ptr v, const _nfloat_fxr_t a, ulong sh, slong w)
{
    slong i, ls = sh / FLINT_BITS;
    unsigned int bs = sh % FLINT_BITS;

    for (i = 0; i < w; i++)
    {
        ulong hi = (i + ls < w) ? a->d[i + ls] : 0;
        ulong hi2 = (i + ls + 1 < w) ? a->d[i + ls + 1] : 0;
        v[i] = (bs == 0) ? hi : ((hi >> bs) | (hi2 << (FLINT_BITS - bs)));
    }
}

/* w = 2 or 3, as constants: the loops unroll, and every limb of v is
   visibly written */
FLINT_FORCE_INLINE void
_fxr_rshift(nn_ptr v, const _nfloat_fxr_t a, ulong sh, slong w)
{
    FLINT_ASSERT(w == 2 || w == 3);

    if (w == 2)
        _fxr_rshift_w(v, a, sh, 2);
    else
        _fxr_rshift_w(v, a, sh, 3);
}

/* r = a + b, positive */
static void
_fxr_add(_nfloat_fxr_t r, const _nfloat_fxr_t a, const _nfloat_fxr_t b, slong w)
{
    ulong u[4], v[4], cy;
    slong e = FLINT_MAX(a->e, b->e);

    _fxr_rshift(u, a, e - a->e, w);
    _fxr_rshift(v, b, e - b->e, w);
    cy = mpn_add_n(u, u, v, w);
    r->err = a->err + b->err + 2;
    if (cy)
    {
        mpn_rshift(u, u, w, 1);
        u[w - 1] |= UWORD(1) << (FLINT_BITS - 1);
        e++;
        r->err += 1;
    }
    flint_mpn_copyi(r->d, u, w);
    r->e = e;
}

/* r = 1/y, y exact */
FLINT_FORCE_INLINE void
_fxr_inv(_nfloat_fxr_t r, const _nfloat_fxr_t y, slong w)
{
    if (w == 2)
    {
        _mp_real_divq_4_2(r->d, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX,
            y->d[1], y->d[0]);
        r->err = 2;
    }
    else
    {
        _mp_real_small_divq_3(&r->d[2], &r->d[1], &r->d[0], (UWORD(1) << (FLINT_BITS - 1)) - 1,
            UWORD_MAX, UWORD_MAX, y->d[2], y->d[1], y->d[0]);
        r->err = 1 + MP_REAL_SMALL_DIVQ_3_ERR;
    }
    /* (1/2 - eps) / D in (1/2, 1) */
    r->e = 1 - y->e;
}

/* the kernels */

FLINT_FORCE_INLINE void
_fx_atan_ratio(nn_ptr r, nn_srcptr y, nn_srcptr x, slong w)
{
    if (w == 2)
        _mp_real_small_atan_ratio_2(&r[1], &r[0], y[1], y[0], x[1], x[0]);
    else
        _mp_real_small_atan_ratio_3(&r[2], &r[1], &r[0], y[2], y[1], y[0], x[2], x[1], x[0]);
}

/* (r, w + 1) = |log(D 2^E)|, returns the sign */
FLINT_FORCE_INLINE int
_fx_log(nn_ptr r, const _nfloat_fxr_t a, slong w)
{
    if (w == 2)
        return _nfloat_small_log_2(&r[2], &r[1], &r[0], a->d[1], a->d[0], a->e);
    else
        return _nfloat_small_log_3(&r[3], &r[2], &r[1], &r[0], a->d[2], a->d[1], a->d[0], a->e);
}

/* (r, w) = log1p(v) for a fraction v < 1 */
FLINT_FORCE_INLINE void
_fx_log1p(nn_ptr r, nn_srcptr v, slong w)
{
    if (w == 2)
        _mp_real_small_log1p_2(&r[1], &r[0], v[1], v[0]);
    else
        _mp_real_small_log1p_3(&r[2], &r[1], &r[0], v[2], v[1], v[0]);
}

/* res = (-1)^sgnbit T 2^k for (T, w + 1) with the integral limb T[w] */
FLINT_FORCE_INLINE int
_fx_set(nfloat_ptr res, nn_srcptr T, slong k, ulong err, int sgnbit, slong w, gr_ctx_t ctx)
{
    if (w == 2)
        return _nfloat_small_set_1(res, T[2], T[1], T[0], k, err, sgnbit, ctx);
    else
        return _nfloat_small_set_2(res, T[3], T[2], T[1], T[0], k, err, sgnbit, ctx);
}

/*
    res = (-1)^sgnbit (c pi/4 + sigma theta), theta = atan2(P, Q) for
    positive P, Q (divided by pi if div_pi), c in {0, 4}, sigma = +-1 with
    sigma = -1 only for c = 4.

    theta = atan(t) for t = P/Q <= 1, else pi/2 - atan(Q/P): with Y the
    smaller and X the larger, t = Y.d 2^-d / X.d (d = X.e - Y.e), within
    2 (Y.err + 1 + X.err) of the exact ratio (|dt| <= (dY + t dX) / X.d,
    X.d >= 1/2), then atan' <= 1. The total c pi/4 + sigma' atan(t) has
    c in {0, 2, 4} (nothing cancels); for c = 0 the relative accuracy
    needs t >= 2^-33 resp. 2^-41.

    Divided by pi: c/4 + sigma' atan(t) (2/pi) / 2, the product truncated
    (below 5 ulps, 2/pi floored: one more) and halved.
*/
static int
_fx_atan2_out(nfloat_ptr res, const _nfloat_fxr_t P, const _nfloat_fxr_t Q,
    int c, int sigma, int div_pi, int sgnbit, slong n, gr_ctx_t ctx)
{
    const _nfloat_fxr_struct * Y;
    const _nfloat_fxr_struct * X;
    ulong y[4], a[4], T[4], err;
    slong w = n + 1, d, i;
    int swap;

    swap = (P->e > Q->e) || (P->e == Q->e && mpn_cmp(P->d, Q->d, w) > 0);
    Y = swap ? Q : P;
    X = swap ? P : Q;
    if (swap)
    {
        c += 2 * sigma;
        sigma = -sigma;
    }

    d = X->e - Y->e;
    if (c == 0 && d > FXR_MAX_SHIFT(n))
        return NFLOAT_FAST_DECLINE;

    err = 2 * (Y->err + 1 + X->err) + ATAN_ERR(w);
    if (d >= FLINT_BITS * w)
    {
        /* t below one unit */
        flint_mpn_zero(a, w);
        err = 2;
    }
    else
    {
        _fxr_rshift(y, Y, d, w);
        _fx_atan_ratio(a, y, X->d, w);
    }

    if (!div_pi)
    {
        if (c == 0)
        {
            a[w] = 0;
            return _fx_set(res, a, 0, err, sgnbit, w, ctx);
        }

        T[w] = mpn_mul_1(T, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
        if (sigma > 0)
            T[w] += mpn_add_n(T, T, a, w);
        else
            T[w] -= mpn_sub_n(T, T, a, w);
        return _fx_set(res, T, 0, err + c, sgnbit, w, ctx);
    }
    else
    {
        nn_srcptr C = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, w);
        ulong p[4];

        if (w == 2)
            _mp_real_mulhi_2x2_sloppy(&p[1], &p[0], a[1], a[0], C[1], C[0]);
        else
            _mp_real_small_mulhi_3x3(&p[2], &p[1], &p[0], a[2], a[1], a[0], C[2], C[1], C[0]);
        mpn_rshift(p, p, w, 1);

        for (i = 0; i <= w; i++)
            T[i] = 0;
        /* c/4 */
        if (c == 4)
            T[w] = 1;
        else
            T[w - 1] = (ulong) c << (FLINT_BITS - 2);
        if (sigma > 0)
            T[w] += mpn_add_n(T, T, p, w);
        else
            T[w] -= mpn_sub_n(T, T, p, w);
        return _fx_set(res, T, 0, err + 8, sgnbit, w, ctx);
    }
}

/* asin, acos (|x| < 1), asec, acsc (|x| > 1) */
static int
_nfloat_asin_acos_fast(nfloat_ptr res, nfloat_srcptr x, int f, int div_pi, gr_ctx_t ctx)
{
    _nfloat_fxr_t Y, U, S, O;
    slong n = NFLOAT_CTX_NLIMBS(ctx), w = n + 1;
    int xsgn = NFLOAT_SGNBIT(x);

    _fxr_set_nfloat(Y, x, n);

    if (f == NFLOAT_ARC_ASIN || f == NFLOAT_ARC_ACOS)
    {
        _fxr_sqr_pm1(U, Y, 0, n, w);
        _fxr_sqrt(S, U, w);
        if (f == NFLOAT_ARC_ASIN)
            return _fx_atan2_out(res, Y, S, 0, 1, div_pi, xsgn, n, ctx);
        else
            return _fx_atan2_out(res, S, Y, xsgn ? 4 : 0, xsgn ? -1 : 1, div_pi, 0, n, ctx);
    }
    else
    {
        _fxr_sqr_pm1(U, Y, -1, n, w);
        _fxr_sqrt(S, U, w);
        _fxr_one(O, w);
        if (f == NFLOAT_ARC_ASEC)
            return _fx_atan2_out(res, S, O, xsgn ? 4 : 0, xsgn ? -1 : 1, div_pi, 0, n, ctx);
        else
            return _fx_atan2_out(res, O, S, 0, 1, div_pi, xsgn, n, ctx);
    }
}

/*
    atanh(r) = (log1p(r) + |log(1 - r)|) / 2 for a fraction (r, w) < 1
    within rerr units: log1p' <= 1, and the factor 1 / (1 - r) <= 2 for
    r <= 1/2 (the callers with rerr > 0); 1 - r is normalized exactly
    (counted as one unit, two in the logarithm), so within
    (LOG1P_ERR + LOG_ERR + 3 rerr + 3) / 2 + 2 units.
*/
static int
_fx_atanh_out(nfloat_ptr res, nn_srcptr r, ulong rerr, int sgnbit, slong w, gr_ctx_t ctx)
{
    _nfloat_fxr_t M;
    ulong u[4], L1[4], L2[4];

    /* 1 - r exactly (as an integer B^w - R, normalized) */
    mpn_neg(u, r, w);
    _fxr_set_mpn(M, u, w, -FLINT_BITS * w, w);

    _fx_log1p(L1, r, w);
    _fx_log(L2, M, w);
    L1[w] = 0;
    mpn_add_n(L2, L2, L1, w + 1);

    return _fx_set(res, L2, -1, (LOG1P_ERR(w) + LOG_ERR(w) + 3 * rerr + 3) / 2 + 2, sgnbit, w, ctx);
}

/*
    asinh(r) for a positive r given as R (exact or within R->err units).

    r >= 1/2: log(A), A = r + sqrt(r^2 + 1) >= 1.6 (r^2 + 1 within one unit
    more than the propagated 2 R->err, the sum within the sum of the errors
    plus 2): the logarithm within LOG_ERR + 2 A.err units (log' = 1/A, A.d
    >= 1/2) of a result >= 0.48.

    r < 1/2: log1p(v), v = r + sqrt(1 + r^2) - 1 < 0.62 as a fraction (from
    S = 2 S.d with S.d at exponent 1), within 2 S.err + R->err + 1 units.

    Only exact R is given with r < 1/2 when r is small (r >= 2^-31 resp.
    2^-39, the result >= 0.9 r).
*/
static int
_fx_asinh_out(nfloat_ptr res, const _nfloat_fxr_t R, int sgnbit, int exact, slong n, gr_ctx_t ctx)
{
    _nfloat_fxr_t U, S, A;
    ulong v[4], s[4], L[4];
    slong w = n + 1;

    if (exact)
        _fxr_sqr_pm1(U, R, 1, n, w);
    else
    {
        /* r^2 + 1 from the w-limb mantissa (n + 1 limbs as "y"), within
           2 R->err + 1 units more */
        ulong z[8], p[8];
        slong zn = 2 * w + 2, K = 2 * FLINT_BITS * w - 2 * R->e;

        if (K < 0 || K >= FLINT_BITS * zn - 1)
            return NFLOAT_FAST_DECLINE;

        flint_mpn_zero(z, zn);
        flint_mpn_zero(p, zn);
        if (w == 2)
            FLINT_MPN_SQR_2X2(z[3], z[2], z[1], z[0], R->d[1], R->d[0]);
        else
            flint_mpn_sqr(z, R->d, 3);
        p[K / FLINT_BITS] = UWORD(1) << (K % FLINT_BITS);
        mpn_add_n(z, z, p, zn);
        _fxr_set_mpn(U, z, zn, 2 * R->e - 2 * FLINT_BITS * w, w);
        U->err += 2 * R->err;
    }

    _fxr_sqrt(S, U, w);

    if (R->e >= 0)
    {
        _fxr_add(A, R, S, w);
        _fx_log(L, A, w);
        return _fx_set(res, L, 0, LOG_ERR(w) + 2 * A->err, sgnbit, w, ctx);
    }
    else
    {
        /* S in [1, 1.12): S.e = 1; S - 1 = 2 S.d - 1 */
        if (S->e != 1)
            return NFLOAT_FAST_DECLINE;
        mpn_lshift(s, S->d, w, 1);
        _fxr_rshift(v, R, -R->e, w);
        mpn_add_n(v, v, s, w);
        _fx_log1p(L, v, w);
        L[w] = 0;
        return _fx_set(res, L, 0, LOG1P_ERR(w) + 2 * S->err + R->err + 2, sgnbit, w, ctx);
    }
}

int
_nfloat_arc_fast(nfloat_ptr res, nfloat_srcptr x, int f, int div_pi, gr_ctx_t ctx)
{
    _nfloat_fxr_t Y, U, S, A, O;
    ulong v[4], L1[4], L2[4];
    slong n = NFLOAT_CTX_NLIMBS(ctx), w = n + 1, e = NFLOAT_EXP(x);
    int xsgn = NFLOAT_SGNBIT(x);

    if (n > 2 || e < FXR_MIN_E(n) || e > FXR_MAX_E)
        return NFLOAT_FAST_DECLINE;

    switch (f)
    {
        case NFLOAT_ARC_ASIN:
        case NFLOAT_ARC_ACOS:
        case NFLOAT_ARC_ASEC:
        case NFLOAT_ARC_ACSC:
            return _nfloat_asin_acos_fast(res, x, f, div_pi, ctx);

        case NFLOAT_ARC_ATAN:
            _fxr_set_nfloat(Y, x, n);
            _fxr_one(O, w);
            return _fx_atan2_out(res, Y, O, 0, 1, div_pi, xsgn, n, ctx);

        case NFLOAT_ARC_ACOT:
            _fxr_set_nfloat(Y, x, n);
            _fxr_one(O, w);
            return _fx_atan2_out(res, O, Y, 0, 1, div_pi, xsgn, n, ctx);

        case NFLOAT_ARC_ASINH:
            _fxr_set_nfloat(Y, x, n);
            return _fx_asinh_out(res, Y, xsgn, 1, n, ctx);

        case NFLOAT_ARC_ACSCH:
            _fxr_set_nfloat(Y, x, n);
            if (e <= 0)
            {
                /* log(1 + S) + |log y|, S = sqrt(1 + y^2) in [1, 1.12) */
                _fxr_sqr_pm1(U, Y, 1, n, w);
                _fxr_sqrt(S, U, w);
                _fxr_one(O, w);
                _fxr_add(A, O, S, w);
                _fx_log(L1, A, w);
                _fx_log(L2, Y, w);
                mpn_add_n(L1, L1, L2, w + 1);
                return _fx_set(res, L1, 0, 2 * LOG_ERR(w) + 2 * A->err, xsgn, w, ctx);
            }
            else
            {
                /* asinh(1/y), 1/y >= 2^-31 resp. 2^-39 */
                if (1 - e < FXR_MIN_E(n))
                    return NFLOAT_FAST_DECLINE;
                _fxr_inv(A, Y, w);
                return _fx_asinh_out(res, A, xsgn, 0, n, ctx);
            }

        case NFLOAT_ARC_ACOSH:
            /* log(A), A = y + sqrt(y^2 - 1): log A >= A - 1 >= S needs
               S >= 2^-31 resp. 2^-39 for the kernel near 1 */
            _fxr_set_nfloat(Y, x, n);
            _fxr_sqr_pm1(U, Y, -1, n, w);
            _fxr_sqrt(S, U, w);
            if (S->e < FXR_MIN_E(n))
                return NFLOAT_FAST_DECLINE;
            _fxr_add(A, Y, S, w);
            _fx_log(L1, A, w);
            return _fx_set(res, L1, 0, LOG_ERR(w) + 2 * A->err, 0, w, ctx);

        case NFLOAT_ARC_ASECH:
            /* log1p(S) + |log y|, S = sqrt(1 - y^2) < 1 as a fraction */
            _fxr_set_nfloat(Y, x, n);
            _fxr_sqr_pm1(U, Y, 0, n, w);
            _fxr_sqrt(S, U, w);
            if (S->e < FXR_MIN_E(n) || S->e > 0)
                return NFLOAT_FAST_DECLINE;
            _fxr_rshift(v, S, -S->e, w);
            _fx_log1p(L1, v, w);
            L1[w] = 0;
            _fx_log(L2, Y, w);
            mpn_add_n(L1, L1, L2, w + 1);
            return _fx_set(res, L1, 0, LOG1P_ERR(w) + LOG_ERR(w) + S->err + 1, 0, w, ctx);

        case NFLOAT_ARC_ATANH:
            /* y < 1 as an exact fraction */
            _fxr_set_nfloat(Y, x, n);
            _fxr_rshift(v, Y, -e, w);
            return _fx_atanh_out(res, v, 0, xsgn, w, ctx);

        default:
            /* acoth, y > 1 */
            _fxr_set_nfloat(Y, x, n);
            if (e == 1)
            {
                /* (log(y + 1) + |log(y - 1)|) / 2, y - 1 and y + 1 exact */
                ulong z[4];

                /* y - 1 = 2 Y.d - 1 as a fraction, nonzero */
                mpn_lshift(z, Y->d, w, 1);
                if (flint_mpn_zero_p(z, w))
                    return NFLOAT_FAST_DECLINE;
                _fxr_set_mpn(U, z, w, -FLINT_BITS * w, w);
                U->err = 0;

                /* y + 1 = (Y.d / 2 + 1/4) 2^2 (exact: the low limb of Y.d
                   is zero) */
                mpn_rshift(z, Y->d, w, 1);
                z[w - 1] += UWORD(1) << (FLINT_BITS - 2);
                flint_mpn_copyi(A->d, z, w);
                A->e = 2;
                A->err = 0;

                _fx_log(L1, A, w);
                _fx_log(L2, U, w);
                mpn_add_n(L1, L1, L2, w + 1);
                return _fx_set(res, L1, -1, (2 * LOG_ERR(w) + 2 * A->err + 1) / 2 + 1, xsgn, w, ctx);
            }
            else
            {
                /* atanh(1/y), 1/y <= 1/2 */
                if (1 - e < FXR_MIN_E(n))
                    return NFLOAT_FAST_DECLINE;
                _fxr_inv(A, Y, w);
                _fxr_rshift(v, A, -A->e, w);
                return _fx_atanh_out(res, v, A->err + 1, xsgn, w, ctx);
            }
    }
}

#endif
