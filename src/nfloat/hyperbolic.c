/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/*
    Hyperbolic functions from E = exp(|x|) and its reciprocal:
    cosh = (E + 1/E) / 2, sinh = (E - 1/E) / 2, tanh = (E - 1/E) / (E + 1/E),
    composed with relative error bounds. For |x| < 1, E - 1/E is computed
    in fixed point with absolute errors at a precision extended by the
    -e + 1 bits cancelling.
*/

/* E = exp(|x|) at w limbs */
static void
_nfloat_exp_abs_rf(_nfloat_rf_t E, nfloat_srcptr x, slong w, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong k;

    k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_EXP(x), 0, w);
    _nfloat_rf_set_mpn(E, y, w + 1, k - FLINT_BITS * w, err, w);
}

/* For 0 < |x| < 1, with E = exp(|x|): sets (T, wb + 1) = E - 1/E and
   (S, wb + 1) = E + 1/E in units of B^-wb, both within the returned
   number of units. The difference is computed in fixed point, so that
   wb limbs with -e extra bits give a relative accuracy for E - 1/E ~ 2x
   (converting E and 1/E to relative errors before subtracting them would
   blow up the error counts). */
static ulong
_nfloat_exp_pm_fixed(nn_ptr T, nn_ptr S, nfloat_srcptr x, slong wb, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2], r[NFLOAT_ELEM_MAX_LIMBS + 2];
    _nfloat_rf_t E, R;
    ulong ey, er;
    slong k;

    k = _nfloat_exp_mpn(y, &ey, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_EXP(x), 0, wb);
    FLINT_ASSERT(k == 0);
    (void) k;

    /* R = 1/E in (1/e, 1), R->exp in {0, -1} */
    _nfloat_rf_set_mpn(E, y, wb + 1, -FLINT_BITS * wb, ey, wb);
    _nfloat_rf_inv(R, E, wb);
    FLINT_ASSERT(R->exp == 0 || R->exp == -1);

    if (R->exp == 0)
    {
        _nfloat_copy_limbs(r, R->d, wb);
        er = R->err;
    }
    else
    {
        mpn_rshift(r, R->d, wb, 1);
        er = R->err + 1;
    }
    r[wb] = 0;

    if (T != NULL)
        mpn_sub_n(T, y, r, wb + 1);
    if (S != NULL)
        mpn_add_n(S, y, r, wb + 1);

    return ey + er;
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs: with |x| = q log 2 + t
    (the reduction of the exp fast paths), exp(|x|) = 2^q Y and
    exp(-|x|) = 2^-q Z, Y = exp(t) and Z = exp(-t) by the paired kernel
    _mp_real_small_exp_pm_2 resp. _3 (one evaluation sharing the
    polynomial, about 1.5 times one exp kernel), so that

        cosh x, |sinh x| = 2^(q-1) (Y +- Z 2^-2q)

    in fixed point. Y is within NFLOAT_SMALL_EXP_ERR units (the kernel,
    the reduction error of t and the dropped limb, times exp' < 2), Z
    <= 1 within the kernel's ERR plus the error of t (below 1.4 units),
    and shifted down by 2q bits within one more unit (or dropped when
    below one unit): 2 ERR + 6 units in all. Without cancellation for
    cosh; sinh keeps 2^-90 resp. 2^-144 relative accuracy for |x| >=
    2^-32 resp. 2^-41 (the callers' ranges).
*/
static int
_nfloat_sinh_cosh_fast_1(nfloat_ptr rs, nfloat_ptr rc, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong xi, x1, x0, q, t2, t1, t0, y2, y1, y0, z2, z1, z0, c2, c1, c0;
    ulong err = 2 * MP_REAL_SMALL_EXP_ERR_2 + 6;
    int status = GR_SUCCESS;

    _nfloat_small_fixed_1(&xi, &x1, &x0, d, e);
    q = _nfloat_small_exp_reduce_3(&t2, &t1, &t0, xi, x1, x0, 0);
    if (2 * q < 2 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_2(&y2, &y1, &y0, &z2, &z1, &z0, t2, t1);
        _nfloat_small_rshift_3(&z2, &z1, &z0, z2, z1, z0, 2 * q);
    }
    else
    {
        /* Z 2^-2q below one unit */
        _mp_real_small_exp_2(&y2, &y1, &y0, t2, t1);
        z2 = z1 = z0 = 0;
    }

    if (rs != NULL)
    {
        sub_dddmmmsss(c2, c1, c0, y2, y1, y0, z2, z1, z0);
        status |= _nfloat_small_set_1(rs, c2, c1, c0, (slong) q - 1, err, sgnbit, ctx);
    }
    if (rc != NULL)
    {
        add_sssaaaaaa(c2, c1, c0, y2, y1, y0, z2, z1, z0);
        status |= _nfloat_small_set_1(rc, c2, c1, c0, (slong) q - 1, err, 0, ctx);
    }
    return status;
}

static int
_nfloat_sinh_cosh_fast_2(nfloat_ptr rs, nfloat_ptr rc, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong xi, x2, x1, x0, q, t3, t2, t1, t0, y3, y2, y1, y0, z3, z2, z1, z0, c3, c2, c1, c0;
    ulong err = 2 * MP_REAL_SMALL_EXP_ERR_3 + 6;
    int status = GR_SUCCESS;

    _nfloat_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);
    q = _nfloat_small_exp_reduce_4(&t3, &t2, &t1, &t0, xi, x2, x1, x0, 0);
    if (2 * q < 3 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_3(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, t3, t2, t1);
        _nfloat_small_rshift_4(&z3, &z2, &z1, &z0, z3, z2, z1, z0, 2 * q);
    }
    else
    {
        _mp_real_small_exp_3(&y3, &y2, &y1, &y0, t3, t2, t1);
        z3 = z2 = z1 = z0 = 0;
    }

    if (rs != NULL)
    {
        sub_ddddmmmmssss(c3, c2, c1, c0, y3, y2, y1, y0, z3, z2, z1, z0);
        status |= _nfloat_small_set_2(rs, c3, c2, c1, c0, (slong) q - 1, err, sgnbit, ctx);
    }
    if (rc != NULL)
    {
        add_ssssaaaaaaaa(c3, c2, c1, c0, y3, y2, y1, y0, z3, z2, z1, z0);
        status |= _nfloat_small_set_2(rc, c3, c2, c1, c0, (slong) q - 1, err, 0, ctx);
    }
    return status;
}

/*
    tanh and coth from the same pair: with N = Y - Z 2^-2q and
    D = Y + Z 2^-2q in fixed point (each within 2 ERR + 6 units, as
    above), |tanh x| = N / D < 1 and |coth x| = D / N > 1 by one division
    in registers. For tanh, N and D are scaled by 2^-sh into [1/2, 1)
    (D in [1, 3); one more unit for each truncation): within
    2 (2 ERR + 6) + 4 + 1 units, the error of D entering at most with
    N / D < 1 and D >= 1. For coth, N = Nn 2^-s with Nn in [1/2, 1) and
    the quotient (D / 8) / Nn < 1/2... at most 3/4: within
    (2 ERR + 6)(2^(max(s, 0) + 1) + 1) + 4 units of the quotient, the
    relative error of N being at most 2^(s+1) (2 ERR + 6) units.
    When Z 2^-2q falls below a unit, the results are 1 within one unit.

    sech and csch (coth = 2, 3) are 2^(1-q) / D resp. 2^(1-q) / N, as
    (1/2 - B^-w) / Dn resp. (1/2 - B^-w) / Nn for the normalized Dn = D 2^-sh
    and Nn = N 2^s: within (2 ERR + 6) + 4 resp.
    (2 ERR + 6) 2^(max(s, 0) + 1) + 4 units of the quotient.
*/
FLINT_FORCE_INLINE ulong
_nfloat_tanh_coth_pair_1(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0,
    ulong d, slong e)
{
    ulong xi, x1, x0, q, t2, t1, t0;

    _nfloat_small_fixed_1(&xi, &x1, &x0, d, e);
    q = _nfloat_small_exp_reduce_3(&t2, &t1, &t0, xi, x1, x0, 0);
    if (2 * q < 2 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_2(y2, y1, y0, z2, z1, z0, t2, t1);
        _nfloat_small_rshift_3(z2, z1, z0, *z2, *z1, *z0, 2 * q);
    }
    else
    {
        _mp_real_small_exp_2(y2, y1, y0, t2, t1);
        *z2 = *z1 = *z0 = 0;
    }
    return q;
}

static int
_nfloat_tanh_coth_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, int coth, gr_ctx_t ctx)
{
    ulong y2, y1, y0, z2, z1, z0, n2, n1, n0, d2, d1, d0, q[2], qq[2], qe;
    ulong e_nd = 2 * MP_REAL_SMALL_EXP_ERR_2 + 6;
    unsigned int sh;

    qe = _nfloat_tanh_coth_pair_1(&y2, &y1, &y0, &z2, &z1, &z0, d, e);

    if (coth <= 1 && z2 == 0 && z1 == 0 && z0 == 0)
        return _nfloat_small_set_1(res, UWORD(1), UWORD(0), UWORD(0), 0, 1, sgnbit, ctx);

    sub_dddmmmsss(n2, n1, n0, y2, y1, y0, z2, z1, z0);
    add_sssaaaaaa(d2, d1, d0, y2, y1, y0, z2, z1, z0);

    if (coth == 2)
    {
        /* sech = 2^(1-q) / D */
        sh = (d2 >= 2) ? 2 : 1;
        _nfloat_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, sh);
        _mp_real_divq_4_2(qq, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX, d1, d0);
        return _nfloat_small_set_1(res, 0, qq[1], qq[0], 2 - (slong) qe - (slong) sh, e_nd + 4, 0, ctx);
    }

    if (coth == 3)
    {
        /* csch = 2^(1-q) / N, N = Nn 2^-s */
        slong s;
        if (n2 != 0)
        {
            n0 = (n0 >> 1) | (n1 << (FLINT_BITS - 1));
            n1 = (n1 >> 1) | (n2 << (FLINT_BITS - 1));
            s = -1;
        }
        else
        {
            sh = flint_clz(n1);
            if (sh != 0)
            {
                n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
                n0 = n0 << sh;
            }
            s = sh;
        }
        _mp_real_divq_4_2(qq, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX, n1, n0);
        return _nfloat_small_set_1(res, 0, qq[1], qq[0], 2 - (slong) qe + s,
            e_nd * (UWORD(2) << FLINT_MAX(s, 0)) + 4, sgnbit, ctx);
    }

    if (!coth)
    {
        /* N / D with both scaled into fractions, D normalized */
        sh = (d2 >= 2) ? 2 : 1;
        _nfloat_small_rshift_3(&n2, &n1, &n0, n2, n1, n0, sh);
        _nfloat_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, sh);
        _mp_real_divq_4_2(q, n1, n0, UWORD(0), UWORD(0), d1, d0);
        return _nfloat_small_set_1(res, 0, q[1], q[0], 0, 2 * e_nd + 5, sgnbit, ctx);
    }
    else
    {
        /* (D / 8) / Nn, N = Nn 2^-s */
        slong s;
        if (n2 != 0)
        {
            /* N in [1, 2) */
            n0 = (n0 >> 1) | (n1 << (FLINT_BITS - 1));
            n1 = (n1 >> 1) | (n2 << (FLINT_BITS - 1));
            s = -1;
        }
        else
        {
            sh = flint_clz(n1);
            if (sh != 0)
            {
                n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
                n0 = n0 << sh;
            }
            s = sh;
        }
        _nfloat_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, 3);
        _mp_real_divq_4_2(q, d1, d0, UWORD(0), UWORD(0), n1, n0);
        return _nfloat_small_set_1(res, 0, q[1], q[0], s + 3,
            e_nd * ((UWORD(2) << FLINT_MAX(s, 0)) + 1) + 4, sgnbit, ctx);
    }
}

FLINT_FORCE_INLINE ulong
_nfloat_tanh_coth_pair_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong d1, ulong d0, slong e)
{
    ulong xi, x2, x1, x0, q, t3, t2, t1, t0;

    _nfloat_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);
    q = _nfloat_small_exp_reduce_4(&t3, &t2, &t1, &t0, xi, x2, x1, x0, 0);
    if (2 * q < 3 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_3(y3, y2, y1, y0, z3, z2, z1, z0, t3, t2, t1);
        _nfloat_small_rshift_4(z3, z2, z1, z0, *z3, *z2, *z1, *z0, 2 * q);
    }
    else
    {
        _mp_real_small_exp_3(y3, y2, y1, y0, t3, t2, t1);
        *z3 = *z2 = *z1 = *z0 = 0;
    }
    return q;
}

static int
_nfloat_tanh_coth_fast_2(nfloat_ptr res, ulong dd1, ulong dd0, slong e, int sgnbit, int coth, gr_ctx_t ctx)
{
    ulong y3, y2, y1, y0, z3, z2, z1, z0, n3, n2, n1, n0, d3, d2, d1, d0, w2, w1, w0, qe;
    ulong e_nd = 2 * MP_REAL_SMALL_EXP_ERR_3 + 6;
    unsigned int sh;

    qe = _nfloat_tanh_coth_pair_2(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, dd1, dd0, e);

    if (coth <= 1 && z3 == 0 && z2 == 0 && z1 == 0 && z0 == 0)
        return _nfloat_small_set_2(res, UWORD(1), UWORD(0), UWORD(0), UWORD(0), 0, 1, sgnbit, ctx);

    sub_ddddmmmmssss(n3, n2, n1, n0, y3, y2, y1, y0, z3, z2, z1, z0);
    add_ssssaaaaaaaa(d3, d2, d1, d0, y3, y2, y1, y0, z3, z2, z1, z0);

    if (coth == 2)
    {
        sh = (d3 >= 2) ? 2 : 1;
        _nfloat_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, sh);
        _mp_real_small_divq_3(&w2, &w1, &w0, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, d2, d1, d0);
        return _nfloat_small_set_2(res, 0, w2, w1, w0, 2 - (slong) qe - (slong) sh,
            e_nd + 3 + MP_REAL_SMALL_DIVQ_3_ERR, 0, ctx);
    }

    if (coth == 3)
    {
        slong s;
        if (n3 != 0)
        {
            _nfloat_small_rshift_4(&n3, &n2, &n1, &n0, n3, n2, n1, n0, 1);
            s = -1;
        }
        else
        {
            sh = flint_clz(n2);
            if (sh != 0)
            {
                n2 = (n2 << sh) | (n1 >> (FLINT_BITS - sh));
                n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
                n0 = n0 << sh;
            }
            s = sh;
        }
        _mp_real_small_divq_3(&w2, &w1, &w0, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, n2, n1, n0);
        return _nfloat_small_set_2(res, 0, w2, w1, w0, 2 - (slong) qe + s,
            e_nd * (UWORD(2) << FLINT_MAX(s, 0)) + 3 + MP_REAL_SMALL_DIVQ_3_ERR, sgnbit, ctx);
    }

    if (!coth)
    {
        sh = (d3 >= 2) ? 2 : 1;
        _nfloat_small_rshift_4(&n3, &n2, &n1, &n0, n3, n2, n1, n0, sh);
        _nfloat_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, sh);
        _mp_real_small_divq_3(&w2, &w1, &w0, n2, n1, n0, d2, d1, d0);
        return _nfloat_small_set_2(res, 0, w2, w1, w0, 0, 2 * e_nd + 4 + MP_REAL_SMALL_DIVQ_3_ERR, sgnbit, ctx);
    }
    else
    {
        slong s;
        if (n3 != 0)
        {
            _nfloat_small_rshift_4(&n3, &n2, &n1, &n0, n3, n2, n1, n0, 1);
            s = -1;
        }
        else
        {
            sh = flint_clz(n2);
            if (sh != 0)
            {
                n2 = (n2 << sh) | (n1 >> (FLINT_BITS - sh));
                n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
                n0 = n0 << sh;
            }
            s = sh;
        }
        _nfloat_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, 3);
        _mp_real_small_divq_3(&w2, &w1, &w0, d2, d1, d0, n2, n1, n0);
        return _nfloat_small_set_2(res, 0, w2, w1, w0, s + 3,
            e_nd * ((UWORD(2) << FLINT_MAX(s, 0)) + 1) + 3 + MP_REAL_SMALL_DIVQ_3_ERR, sgnbit, ctx);
    }
}

/* tanh (coth = 0), coth (coth = 1), sech (2) or csch (3) of a normal x by
   the fast paths: for
   2^-32 <= |x| < 2^61 (resp. 2^-41), NFLOAT_FAST_DECLINE otherwise */
int
_nfloat_tanh_coth_fast(nfloat_ptr res, nfloat_srcptr x, int coth, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), e = NFLOAT_EXP(x);

    if (n > 2 || e > NFLOAT_EXP_MAX_ARG_EXP || e < ((n == 1) ? -31 : -40))
        return NFLOAT_FAST_DECLINE;

    if (n == 1)
        return _nfloat_tanh_coth_fast_1(res, NFLOAT_D(x)[0], e, NFLOAT_SGNBIT(x), coth, ctx);
    else
        return _nfloat_tanh_coth_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, NFLOAT_SGNBIT(x), coth, ctx);
}
#endif

/* sinh (want & 1), cosh (want & 2) */
static int
_nfloat_sinh_cosh(nfloat_ptr rs, nfloat_ptr rc, nfloat_srcptr x, gr_ctx_t ctx)
{
    _nfloat_rf_t E, R, T;
    slong n, e, w;
    int sgnbit, status = GR_SUCCESS;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
    {
        if (rs != NULL)
            status |= _nfloat_overflow(rs, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_overflow(rc, 0, ctx);
        return status;
    }

    /* sinh x = x (1 + delta), cosh x = 1 + delta', delta, delta' < x^2 */
    if (2 * e <= -FLINT_BITS * n)
    {
        if (rs != NULL)
            status |= _nfloat_set_tiny_rel_err(rs, x, 1, ctx);
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 0, -FLINT_BITS * n, ctx);
        return status;
    }

#if FLINT_BITS == 64
    if (n == 1 && (rs == NULL || e >= -31))
        return _nfloat_sinh_cosh_fast_1(rs, rc, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2 && (rs == NULL || e >= -40))
        return _nfloat_sinh_cosh_fast_2(rs, rc, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    if (e <= 0)
    {
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 2], S[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong err;

        w = _nfloat_elem_limbs(ctx, -e + 1);
        err = _nfloat_exp_pm_fixed((rs != NULL) ? T : NULL, (rc != NULL) ? S : NULL, x, w, ctx);

        if (rs != NULL)
            status |= _nfloat_set_mpn_err(rs, T, w + 1, -FLINT_BITS * w - 1, err, err, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_set_mpn_err(rc, S, w + 1, -FLINT_BITS * w - 1, err, err, 0, ctx);
        return status;
    }

    w = _nfloat_elem_limbs(ctx, 0);

    _nfloat_exp_abs_rf(E, x, w, ctx);
    _nfloat_rf_inv(R, E, w);

    if (rs != NULL)
    {
        _nfloat_rf_sub(T, E, R, w);
        T->exp--;
        status |= _nfloat_set_rf(rs, T, sgnbit, w, ctx);
    }

    if (rc != NULL)
    {
        _nfloat_rf_add(T, E, R, w);
        T->exp--;
        status |= _nfloat_set_rf(rc, T, 0, w, ctx);
    }

    return status;
}

int
nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_set(res, x, ctx) | (NFLOAT_IS_NAN(x) ? nfloat_nan(res, ctx) : GR_SUCCESS);

    return _nfloat_sinh_cosh(res, NULL, x, ctx);
}

int
nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_sinh_cosh(NULL, res, x, ctx);
}

int
nfloat_sinh_cosh(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        /* sinh first (it leaves special values unchanged, so res1 may
           alias x), as two statements: the operands of | are
           unsequenced */
        int status = nfloat_sinh(res1, x, ctx);
        return status | nfloat_cosh(res2, x, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_sinh_cosh(res1, res2, t, ctx);
    }

    return _nfloat_sinh_cosh(res1, res2, x, ctx);
}

int
nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    _nfloat_rf_t E, R, N, D;
    slong n, e, w;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* tanh x = x (1 - delta), delta < x^2 / 3 */
    if (2 * e <= -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, -1, ctx);

#if FLINT_BITS == 64
    {
        int status = _nfloat_tanh_coth_fast(res, x, 0, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif

    /* |tanh x| = 1 - 2 / (E^2 + 1) with 2 / (E^2 + 1) < 2 exp(-2 |x|) <
       2^(1 - 2.88 |x|) <= 2^(-FLINT_BITS n) for |x| >= 2^(e - 1) >= 23 n + 1 */
    if (e - 1 >= (slong) FLINT_BIT_COUNT(23 * n + 1))
        return _nfloat_one_plus_tiny(res, sgnbit, 1, -FLINT_BITS * n, ctx);

    w = _nfloat_elem_limbs(ctx, 0);

    if (e <= 0)
    {
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 2], S[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong err;
        /* E - 1/E in [2^e, 2^(e+2)) B^wb, with -e + 1 more bits for the
           cancellation, and enough that its normalization to w limbs
           shifts left by at most NFLOAT_ELEM_MAX_NORM_SHIFT bits */
        slong wb = _nfloat_elem_limbs(ctx, -e + 1);
        if (-e > NFLOAT_ELEM_MAX_NORM_SHIFT)
            wb = FLINT_MAX(wb, w + (-e - NFLOAT_ELEM_MAX_NORM_SHIFT + FLINT_BITS - 1) / FLINT_BITS);

        /* (E - 1/E) / (E + 1/E), the difference in fixed point */
        err = _nfloat_exp_pm_fixed(T, S, x, wb, ctx);
        _nfloat_rf_set_mpn(N, T, wb + 1, -FLINT_BITS * wb, err, w);
        _nfloat_rf_set_mpn(D, S, wb + 1, -FLINT_BITS * wb, err, w);
    }
    else
    {
        _nfloat_exp_abs_rf(E, x, w, ctx);
        _nfloat_rf_inv(R, E, w);
        _nfloat_rf_sub(N, E, R, w);
        _nfloat_rf_add(D, E, R, w);
    }

    _nfloat_rf_div(N, N, D, w);
    return _nfloat_set_rf(res, N, sgnbit, w, ctx);
}
