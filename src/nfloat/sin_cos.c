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
    Trigonometric functions.

    For |x| >= 1, |x| is reduced modulo pi/2 in fixed point by the
    reduction of mp_real_sin_cos_bits (_mp_real_trig_reduce, on the
    mantissa realigned to the limb radix): |x| = a pi/2 + s v, s = +-1,
    v in [0, pi/4] at w fraction limbs within 3 ulps, and

        sin |x| = s sin v, cos v, -s sin v, -cos v   for a = 0, 1, 2, 3,
        cos |x| = cos v, -s sin v, -cos v, s sin v,
        tan |x| = s tan v (a even), -s / tan v (a odd).

    When the needed output is sin v or tan v and v has z leading zero
    bits, the reduction is repeated with z more bits so that the result
    keeps its relative accuracy. Arguments below 1 go to the kernels
    unreduced.

    sin(pi x), cos(pi x), tan(pi x) reduce |x| exactly by the reduction of
    mp_real_sin_cos_pi_bits (_mp_real_trig_pi_reduce): |x| = a/2 + s t
    with t in [0, 1/4], then v = pi t at enough limbs for t's leading
    zero bits.
*/

/* native reduction up to this exponent; beyond, mp_real balls */
#define NFLOAT_TRIG_MAX_EXP 65536

/* sin, cos, tan (op = 0, 1, 2) and of pi x (op = 3, 4, 5) by mp_real
   balls, for huge arguments and for more cancellation than the
   temporaries allow */
static int
_nfloat_trig_ball(mp_real_t r, slong * shift, const mp_real_t x, const mp_real_t y, int op, slong prec)
{
    (void) y;
    *shift = 0;

    switch (op)
    {
        case 0: mp_real_sin_cos_bits(r, NULL, x, prec); return 1;
        case 1: mp_real_sin_cos_bits(NULL, r, x, prec); return 1;
        case 2: return mp_real_tan_bits(r, x, prec);
        case 3: mp_real_sin_cos_pi_bits(r, NULL, x, prec); return 1;
        case 4: mp_real_sin_cos_pi_bits(NULL, r, x, prec); return 1;
        default: return mp_real_tan_pi_bits(r, x, prec);
    }
}

/* the outputs that are not NULL by the above (pi: 0 or 3) */
static int
_nfloat_trig_slow(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, int pi, gr_ctx_t ctx)
{
    slong min_prec = pi ? 0 : _nfloat_trig_ball_min_prec(x);
    int status = GR_SUCCESS;

    if (rs != NULL)
        status |= _nfloat_eval_ball(rs, _nfloat_trig_ball, x, NULL, pi + 0, min_prec, ctx);
    if (rc != NULL)
        status |= _nfloat_eval_ball(rc, _nfloat_trig_ball, x, NULL, pi + 1, min_prec, ctx);
    if (rt != NULL)
        status |= _nfloat_eval_ball(rt, _nfloat_trig_ball, x, NULL, pi + 2, min_prec, ctx);

    return status;
}

/* the sine and cosine results of the kernels (ys, yc at w fraction
   limbs within err ulps, with the signs code) into rs, rc (either may
   be NULL); the sine of x is odd */
static int
_nfloat_trig_set(nfloat_ptr rs, nfloat_ptr rc, nn_ptr ys, nn_ptr yc, slong w, ulong err, int code, int xsgn, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (rs != NULL)
        status |= _nfloat_set_mpn_err(rs, ys, w + 1, -FLINT_BITS * w, err, err, (code & 1) ^ xsgn, ctx);
    if (rc != NULL)
        status |= _nfloat_set_mpn_err(rc, yc, w + 1, -FLINT_BITS * w, err, err, (code >> 1) & 1, ctx);

    return status;
}

#if FLINT_BITS == 64
/*
    In-register fast paths for sin, cos and tan and for their versions
    at pi x, at one and two limbs (64-bit only). The reductions give
    |x| = a pi/2 + s v resp. pi |x| = a pi/2 + s v with v in [0, pi/4]
    as a fixed-point number for the table-driven sine and cosine kernel
    of mp_real at two resp. three limbs, within ev ulps of the kernel:

    * x itself for 2^-33 <= |x| < 2^63 (resp. 2^-41): |x| < 1 goes to the
      kernel unreduced (exactly, ev = 0) except for the tangent of
      |x| >= 1/2, larger |x| is reduced with q
      from a one-limb 2/pi (one low at most, corrected once) and 2q times
      the floor P of pi/4 at three resp. four fraction limbs (an error of
      2(q + 1) units of B^-3 resp. B^-4, below 0.65 ulp of B^-2 resp.
      B^-3 for q < 2^62.4), folded to v in [0, pi/4]; the kernel input
      drops the last limb (ev = 2).

    * pi x for 2^-32 <= |x| < 2^63 (resp. 2^-40): 2|x| = m + f exactly,
      t = f/2 or (1 - f)/2 in [0, 1/4] (exact at three limbs, or within
      half an ulp at two limbs), and v = 3t + (pi - 3) t (the product
      truncated, pi - 3 floored: within 6 ulps of B^-3; ev = 2 resp. 8).
      Half-integers are left to the general code (exact results).

    The sine and cosine of v give sin and cos of |x| by the assignment of
    _nfloat_trig_eval; the tangent and the cotangent take tan v = sin v /
    cos v (one division: 3 (ERR + ev) + 1 units with sin v < 1/sqrt(2) <
    cos v) or cot v = cos v / sin v (as cos v / (2 S) with S = sin v 2^sh
    normalized: (ERR + ev)(2^(sh+1) + 2) + 1 units, sin v >= 2^-(sh+1)),
    for a even resp. odd (tan) or odd resp. even (cot). The secant and
    the cosecant take 1 / D for D = cos v or sin v (by the same
    assignment), as (1/2 - B^-4) / (D 2^sh) with D 2^sh normalized:
    (ERR + ev) 2^(sh+1) + 2 units. sinc divides V 2^sh / 2 (V = sin |x|
    by the same assignment) by the mantissa X of x, or by X pi/4
    normalized for the function of pi x: (ERR + ev) 2^(sh+1) + 2 units,
    plus 16 (one limb) resp. 24 (two limbs) for the relative error of
    X pi/4.
    At three limbs the division (_mp_real_small_divq_3) adds
    MP_REAL_SMALL_DIVQ_3_ERR units instead of one.

    The fast paths decline (return NFLOAT_TRIG_FAST_DECLINE) when v <
    2^-32 resp. 2^-40, where sin v (or tan v) would lose too much
    relative accuracy, and for tan v with sin v >= cos v (v within a few
    ulps of pi/4).
*/

#define NFLOAT_TRIG_FAST_DECLINE NFLOAT_FAST_DECLINE

/* the outputs from sin v, cos v as fixed-point numbers Y 2^(-FLINT_BITS w)
   by the assignment of _nfloat_trig_eval */
/* the outputs from v at two limbs (v within ev ulps) */
static int
_nfloat_trig_out_fast_1(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, int cot, int a, int sneg,
    ulong v1, ulong v0, ulong ev, int xsgn, nfloat_srcptr x, int pi, gr_ctx_t ctx)
{
    ulong s1, s0, g1, g0, c2, c1, c0, q[2];
    ulong err = MP_REAL_SMALL_SIN_COS_ERR_2 + ev;
    int nsv, ncv, status = GR_SUCCESS;
    nfloat_ptr os, oc;

    _mp_real_small_sin_cos_g_2(&s1, &s0, &g1, &g0, v1, v0);

    if (rt != NULL)
    {
        /* cos v = 1 - g < 1 (g > 0 for v >= 2^-32) */
        sub_ddmmss(c1, c0, UWORD(0), UWORD(0), g1, g0);

        if (cot == NFLOAT_TRIG_SINC)
        {
            /* sin |x| / |x| (or / (pi |x|)): V = sin |x| from s (a even)
               or c (a odd), V 2^sh normalized, halved, divided by
               X = the mantissa of x (or by P = X pi/4 2^psh normalized) */
            ulong V1 = (a & 1) ? c1 : s1, V0 = (a & 1) ? c0 : s0, X1, X0;
            unsigned int sh = flint_clz(V1), psh = 0;
            ulong perr = 0;
            int sg;

            ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
            sg = (a & 1) ? ncv : nsv;
            if (sh != 0)
            {
                V1 = (V1 << sh) | (V0 >> (FLINT_BITS - sh));
                V0 = V0 << sh;
            }
            X1 = NFLOAT_D(x)[0];
            X0 = 0;
            if (pi)
            {
                nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 2);
                _mp_real_mulhi_2x2_sloppy(&X1, &X0, X1, X0, P[1], P[0]);
                psh = !(X1 >> (FLINT_BITS - 1));
                if (psh)
                {
                    X1 = (X1 << 1) | (X0 >> (FLINT_BITS - 1));
                    X0 = X0 << 1;
                }
                /* within 4 units, 8 normalized: relative 16 */
                perr = 16;
            }
            _mp_real_divq_4_2(q, V1 >> 1, (V1 << (FLINT_BITS - 1)) | (V0 >> 1), V0 << (FLINT_BITS - 1), UWORD(0), X1, X0);
            return _nfloat_small_set_1(rt, 0, q[1], q[0],
                (pi ? -1 + (slong) psh : 1) - (slong) sh - NFLOAT_EXP(x),
                err * (UWORD(2) << sh) + perr + 2, sg, ctx);
        }

        if (cot >= NFLOAT_TRIG_SEC)
        {
            /* sec |x| = +-1 / cos v (a even) or +-1 / sin v (a odd), csc |x|
               conversely; D = the denominator 2^sh, normalized */
            int dis = (a & 1) ^ (cot == NFLOAT_TRIG_CSC);
            ulong D1 = dis ? s1 : c1, D0 = dis ? s0 : c0;
            unsigned int sh = flint_clz(D1);
            int sg;

            ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
            sg = (dis ? nsv : ncv) ^ ((cot == NFLOAT_TRIG_CSC) ? xsgn : 0);
            if (sh != 0)
            {
                D1 = (D1 << sh) | (D0 >> (FLINT_BITS - sh));
                D0 = D0 << sh;
            }
            _mp_real_divq_4_2(q, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX, D1, D0);
            return _nfloat_small_set_1(rt, 0, q[1], q[0], sh + 1, err * (UWORD(2) << sh) + 2, sg, ctx);
        }

        /* tan |x| = s tan v (a even), -s cot v (a odd); cot |x| = s cot v
           (a even), -s tan v (a odd) */
        if (!((a & 1) ^ cot))
        {
            /* tan v = sin v / cos v < 1 */
            if (s1 > c1 || (s1 == c1 && s0 >= c0))
                return NFLOAT_TRIG_FAST_DECLINE;
            _mp_real_divq_4_2(q, s1, s0, UWORD(0), UWORD(0), c1, c0);
            return _nfloat_small_set_1(rt, 0, q[1], q[0], 0, 3 * err + 1, (a & 1) ^ sneg ^ xsgn, ctx);
        }
        else
        {
            /* cot v = 2^(sh+1) cos v / (2 S), S = sin v 2^sh */
            unsigned int sh = flint_clz(s1);
            if (sh != 0)
            {
                s1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                s0 = s0 << sh;
            }
            _mp_real_divq_4_2(q, c1 >> 1, (c1 << (FLINT_BITS - 1)) | (c0 >> 1),
                c0 << (FLINT_BITS - 1), UWORD(0), s1, s0);
            return _nfloat_small_set_1(rt, 0, q[1], q[0], sh + 1,
                err * ((UWORD(2) << sh) + 2) + 1, (a & 1) ^ sneg ^ xsgn, ctx);
        }
    }

    ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
    os = (a & 1) ? rc : rs;
    oc = (a & 1) ? rs : rc;

    if (os != NULL)
        status |= _nfloat_small_set_1(os, 0, s1, s0, 0, err,
            nsv ^ ((os == rs) ? xsgn : 0), ctx);
    if (oc != NULL)
    {
        sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
        status |= _nfloat_small_set_1(oc, c2, c1, c0, 0, err,
            ncv ^ ((oc == rs) ? xsgn : 0), ctx);
    }

    return status;
}

static int
_nfloat_trig_out_fast_2(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, int cot, int a, int sneg,
    ulong v2, ulong v1, ulong v0, ulong ev, int xsgn, nfloat_srcptr x, int pi, gr_ctx_t ctx)
{
    ulong s2, s1, s0, g2, g1, g0, c3, c2, c1, c0, w2, w1, w0;
    ulong err = MP_REAL_SMALL_SIN_COS_ERR_3 + ev;
    int nsv, ncv, status = GR_SUCCESS;
    nfloat_ptr os, oc;

    _mp_real_small_sin_cos_g_3(&s2, &s1, &s0, &g2, &g1, &g0, v2, v1, v0);

    if (rt != NULL)
    {
        sub_dddmmmsss(c2, c1, c0, UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);

        if (cot == NFLOAT_TRIG_SINC)
        {
            ulong V2 = (a & 1) ? c2 : s2, V1 = (a & 1) ? c1 : s1, V0 = (a & 1) ? c0 : s0;
            ulong X2, X1, X0;
            unsigned int sh = flint_clz(V2), psh = 0;
            ulong perr = 0;
            int sg;

            ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
            sg = (a & 1) ? ncv : nsv;
            if (sh != 0)
            {
                V2 = (V2 << sh) | (V1 >> (FLINT_BITS - sh));
                V1 = (V1 << sh) | (V0 >> (FLINT_BITS - sh));
                V0 = V0 << sh;
            }
            /* halved, the lowest bit dropped (half a unit) */
            V0 = (V0 >> 1) | (V1 << (FLINT_BITS - 1));
            V1 = (V1 >> 1) | (V2 << (FLINT_BITS - 1));
            V2 = V2 >> 1;
            X2 = NFLOAT_D(x)[1];
            X1 = NFLOAT_D(x)[0];
            X0 = 0;
            if (pi)
            {
                nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
                _mp_real_small_mulhi_3x3(&X2, &X1, &X0, X2, X1, X0, P[2], P[1], P[0]);
                psh = !(X2 >> (FLINT_BITS - 1));
                if (psh)
                {
                    X2 = (X2 << 1) | (X1 >> (FLINT_BITS - 1));
                    X1 = (X1 << 1) | (X0 >> (FLINT_BITS - 1));
                    X0 = X0 << 1;
                }
                /* within 6 units, 12 normalized: relative 24 */
                perr = 24;
            }
            _mp_real_small_divq_3(&w2, &w1, &w0, V2, V1, V0, X2, X1, X0);
            return _nfloat_small_set_2(rt, 0, w2, w1, w0,
                (pi ? -1 + (slong) psh : 1) - (slong) sh - NFLOAT_EXP(x),
                err * (UWORD(2) << sh) + perr + 2 + MP_REAL_SMALL_DIVQ_3_ERR, sg, ctx);
        }

        if (cot >= NFLOAT_TRIG_SEC)
        {
            int dis = (a & 1) ^ (cot == NFLOAT_TRIG_CSC);
            ulong D2 = dis ? s2 : c2, D1 = dis ? s1 : c1, D0 = dis ? s0 : c0;
            unsigned int sh = flint_clz(D2);
            int sg;

            ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
            sg = (dis ? nsv : ncv) ^ ((cot == NFLOAT_TRIG_CSC) ? xsgn : 0);
            if (sh != 0)
            {
                D2 = (D2 << sh) | (D1 >> (FLINT_BITS - sh));
                D1 = (D1 << sh) | (D0 >> (FLINT_BITS - sh));
                D0 = D0 << sh;
            }
            _mp_real_small_divq_3(&w2, &w1, &w0, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, D2, D1, D0);
            return _nfloat_small_set_2(rt, 0, w2, w1, w0, sh + 1,
                err * (UWORD(2) << sh) + 1 + MP_REAL_SMALL_DIVQ_3_ERR, sg, ctx);
        }

        if (!((a & 1) ^ cot))
        {
            if (s2 > c2 || (s2 == c2 && (s1 > c1 || (s1 == c1 && s0 >= c0))))
                return NFLOAT_TRIG_FAST_DECLINE;
            _mp_real_small_divq_3(&w2, &w1, &w0, s2, s1, s0, c2, c1, c0);
            return _nfloat_small_set_2(rt, 0, w2, w1, w0, 0,
                3 * err + MP_REAL_SMALL_DIVQ_3_ERR, (a & 1) ^ sneg ^ xsgn, ctx);
        }
        else
        {
            unsigned int sh = flint_clz(s2);
            if (sh != 0)
            {
                s2 = (s2 << sh) | (s1 >> (FLINT_BITS - sh));
                s1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                s0 = s0 << sh;
            }
            /* cos v / 2, the lowest bit dropped (half a unit) */
            _mp_real_small_divq_3(&w2, &w1, &w0, c2 >> 1, (c2 << (FLINT_BITS - 1)) | (c1 >> 1),
                (c1 << (FLINT_BITS - 1)) | (c0 >> 1), s2, s1, s0);
            return _nfloat_small_set_2(rt, 0, w2, w1, w0, sh + 1,
                err * ((UWORD(2) << sh) + 2) + MP_REAL_SMALL_DIVQ_3_ERR + 1, (a & 1) ^ sneg ^ xsgn, ctx);
        }
    }

    ARF_SMALL_TRIG_SIGNS(nsv, ncv, a, sneg, xsgn);
    os = (a & 1) ? rc : rs;
    oc = (a & 1) ? rs : rc;

    if (os != NULL)
        status |= _nfloat_small_set_2(os, 0, s2, s1, s0, 0, err,
            nsv ^ ((os == rs) ? xsgn : 0), ctx);
    if (oc != NULL)
    {
        sub_ddddmmmmssss(c3, c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
        status |= _nfloat_small_set_2(oc, c3, c2, c1, c0, 0, err,
            ncv ^ ((oc == rs) ? xsgn : 0), ctx);
    }

    return status;
}

/* sin, cos, tan of x (pi = 0) or of pi x (pi = 1) */
static int
_nfloat_trig_fast(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, int pi, int cot, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), e = NFLOAT_EXP(x);
    int xsgn = NFLOAT_SGNBIT(x), a = 0, sneg = 0;
    ulong v2, v1, v0, ev;

    if (e >= FLINT_BITS)
        return NFLOAT_TRIG_FAST_DECLINE;

    if (n == 1)
    {
        ulong d = NFLOAT_D(x)[0];

        if (pi)
        {
            if (e < -31 || _arf_small_trig_pi_reduce_1(&a, &sneg, &v1, &v0, d, e) != 0)
                return NFLOAT_TRIG_FAST_DECLINE;
            ev = 2;
        }
        else if (e < 0 || (e == 0 && rt == NULL))
        {
            /* (the tangent folds |x| > pi/4) */
            if (e < -32)
                return NFLOAT_TRIG_FAST_DECLINE;
            v1 = (e == 0) ? d : (d >> (-e));
            v0 = (e == 0) ? 0 : (d << (FLINT_BITS + e));
            ev = 0;
        }
        else
        {
            if (_arf_small_trig_reduce_1(&a, &sneg, &v1, &v0, d, e) != 0)
                return NFLOAT_TRIG_FAST_DECLINE;
            ev = 2;
        }

        return _nfloat_trig_out_fast_1(rs, rc, rt, cot, a, sneg, v1, v0, ev, xsgn, x, pi, ctx);
    }
    else
    {
        ulong d1 = NFLOAT_D(x)[1], d0 = NFLOAT_D(x)[0];

        if (pi)
        {
            if (e < -39 || _arf_small_trig_pi_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
                return NFLOAT_TRIG_FAST_DECLINE;
            ev = 8;
        }
        else if (e < 0 || (e == 0 && rt == NULL))
        {
            if (e < -40)
                return NFLOAT_TRIG_FAST_DECLINE;
            if (e == 0)
            {
                v2 = d1; v1 = d0; v0 = 0;
            }
            else
            {
                unsigned int sh = -e;
                v2 = d1 >> sh;
                v1 = (d1 << (FLINT_BITS - sh)) | (d0 >> sh);
                v0 = d0 << (FLINT_BITS - sh);
            }
            ev = 0;
        }
        else
        {
            if (_arf_small_trig_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
                return NFLOAT_TRIG_FAST_DECLINE;
            ev = 2;
        }

        return _nfloat_trig_out_fast_2(rs, rc, rt, cot, a, sneg, v2, v1, v0, ev, xsgn, x, pi, ctx);
    }
}
#endif

#if FLINT_BITS == 64
/* cot, sec, csc or sinc (which = NFLOAT_TRIG_COT, _SEC, _CSC, _SINC) of x
   (pi = 0) or of pi x (pi = 1) by the fast paths; NFLOAT_FAST_DECLINE if
   they do not apply */
int
_nfloat_trig_recip_fast(nfloat_ptr res, nfloat_srcptr x, int pi, int which, gr_ctx_t ctx)
{
    if (NFLOAT_CTX_NLIMBS(ctx) > 2)
        return NFLOAT_FAST_DECLINE;
    return _nfloat_trig_fast(NULL, NULL, res, x, pi, which, ctx);
}
#endif

/* which = 0: sin and/or cos; which = 1: tan */
static int
_nfloat_trig(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong ys[NFLOAT_ELEM_MAX_LIMBS + 2], yc[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, w, pb;
    int xsgn;
    int status;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    xsgn = NFLOAT_SGNBIT(x);

    /* tiny x: sin x = x (1 - x^2/6 + ...), cos x = 1 - x^2/2 + ...,
       tan x = x (1 + x^2/3 + ...), x^2 < 2^(2e) <= 2^(-FLINT_BITS n) */
    if (2 * e <= -FLINT_BITS * n)
    {
        status = GR_SUCCESS;
        if (rt != NULL)
            return _nfloat_set_tiny_rel_err(rt, x, 1, ctx);
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 1, -FLINT_BITS * n, ctx);
        if (rs != NULL)
            status |= _nfloat_set_tiny_rel_err(rs, x, -1, ctx);
        return status;
    }

#if FLINT_BITS == 64
    if (n <= 2)
    {
        status = _nfloat_trig_fast(rs, rc, rt, x, 0, 0, ctx);
        if (status != NFLOAT_TRIG_FAST_DECLINE)
            return status;
    }
#endif

    if (e > NFLOAT_TRIG_MAX_EXP)
        return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);

    pb = _nfloat_elem_bits(ctx, 0);

    if (rt == NULL)
    {
        int code = _mp_real_sin_cos_fp(ys, yc, &w, &err, NFLOAT_D(x), n, e, pb, NFLOAT_ELEM_MAX_LIMBS, (rs != NULL) | ((rc != NULL) << 1));

        /* more cancellation than the temporaries allow */
        if (w == NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);

        return _nfloat_trig_set(rs, rc, ys, yc, w, err, code, xsgn, ctx);
    }
    else
    {
        int recip, neg = _mp_real_tan_fp(ys, &w, &err, &recip, NFLOAT_D(x), n, e, pb, NFLOAT_ELEM_MAX_LIMBS);

        if (w == NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);

        return _nfloat_set_mpn_err_recip(rt, ys, w + 1, -FLINT_BITS * w, err, recip, neg ^ xsgn, ctx);
    }
}

int
nfloat_sin_cos(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res1, ctx) | nfloat_one(res2, ctx);
        return nfloat_nan(res1, ctx) | nfloat_nan(res2, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_trig(res1, res2, NULL, t, ctx);
    }

    return _nfloat_trig(res1, res2, NULL, x, ctx);
}

int
nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(res, NULL, NULL, x, ctx);
}

int
nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(NULL, res, NULL, x, ctx);
}

int
nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(NULL, NULL, res, x, ctx);
}

/* sin(pi x), cos(pi x), tan(pi x) *******************************************/

static int
_nfloat_trig_pi(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong ys[2 * NFLOAT_ELEM_MAX_LIMBS + 4], yc[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, w, k, pb;
    int xsgn;
    int status = GR_SUCCESS;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    xsgn = NFLOAT_SGNBIT(x);

    /* an even integer (all bits of weight >= 2) */
    if (e > FLINT_BITS * n)
    {
        if (rs != NULL)
            status |= nfloat_zero(rs, ctx);
        if (rc != NULL)
            status |= nfloat_one(rc, ctx);
        if (rt != NULL)
            status |= nfloat_zero(rt, ctx);
        return status;
    }

    /* tiny x: sin(pi x) = pi x (1 - delta), tan(pi x) = pi x (1 + delta),
       cos(pi x) = 1 - delta', delta < (pi x)^2 < 2^(2e + 4) */
    if (2 * e + 4 <= -FLINT_BITS * n - 2 * FLINT_BITS)
    {
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 1, -FLINT_BITS * n, ctx);

        if (rs != NULL || rt != NULL)
        {
            /* pi |x| from (d, n) times pi at n + 2 limbs: within 3 units of
               B^-(n + 1) relative to the mantissa, delta below one more */
            ulong t[2 * NFLOAT_MAX_LIMBS + 4];
            ulong p[NFLOAT_MAX_LIMBS + 3];
            slong pn = n + 2;

            p[pn] = mpn_lshift(p, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, pn), pn, 2);
            flint_mpn_mul(t, p, pn + 1, NFLOAT_D(x), n);
            /* (t, 2n + 3) = pi D B^pn (to 4 ulps of the pi limbs); keep the
               top n + 2 limbs: units of B^(n + 1) */
            if (rs != NULL)
                status |= _nfloat_set_mpn_err(rs, t + n + 1, n + 2,
                    e - FLINT_BITS * n + FLINT_BITS * (n + 1) - FLINT_BITS * pn, 6, 6, xsgn, ctx);
            if (rt != NULL)
                status |= _nfloat_set_mpn_err(rt, t + n + 1, n + 2,
                    e - FLINT_BITS * n + FLINT_BITS * (n + 1) - FLINT_BITS * pn, 6, 6, xsgn, ctx);
        }
        return status;
    }

#if FLINT_BITS == 64
    if (n <= 2)
    {
        status = _nfloat_trig_fast(rs, rc, rt, x, 1, 0, ctx);
        if (status != NFLOAT_TRIG_FAST_DECLINE)
            return status;
        status = GR_SUCCESS;
    }
#endif

    pb = _nfloat_elem_bits(ctx, 0);

    if (rt == NULL)
    {
        int code = _mp_real_sin_cos_pi_fp(ys, yc, &w, &k, &err, NFLOAT_D(x), n, e, pb, NFLOAT_ELEM_MAX_LIMBS, (rs != NULL) | ((rc != NULL) << 1));

        FLINT_ASSERT(!(code & MP_REAL_TRIG_TINY));

        if (code & MP_REAL_TRIG_EXACT)
        {
            /* x a multiple of 1/2: sin, cos in {0, 1, -1} */
            int a = code & 3;

            if (a & 1)
            {
                if (rs != NULL)
                {
                    status |= nfloat_one(rs, ctx);
                    NFLOAT_SGNBIT(rs) = (a == 3) ^ xsgn;
                }
                if (rc != NULL)
                    status |= nfloat_zero(rc, ctx);
            }
            else
            {
                if (rs != NULL)
                    status |= nfloat_zero(rs, ctx);
                if (rc != NULL)
                {
                    status |= nfloat_one(rc, ctx);
                    NFLOAT_SGNBIT(rc) = (a == 2);
                }
            }
            return status;
        }

        if (w == NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 3, ctx);

        return _nfloat_trig_set(rs, rc, ys, yc, w, err, code, xsgn, ctx);
    }
    else
    {
        int recip, code = _mp_real_tan_pi_fp(ys, &w, &k, &err, &recip, NFLOAT_D(x), n, e, pb, NFLOAT_ELEM_MAX_LIMBS);

        FLINT_ASSERT(!(code & MP_REAL_TRIG_TINY));

        if (code & MP_REAL_TRIG_EXACT)
        {
            if (code & MP_REAL_TRIG_ONE)
            {
                status = nfloat_one(rt, ctx);
                NFLOAT_SGNBIT(rt) = (code & 1) ^ xsgn;
                return status;
            }
            /* a multiple of 1/2: 0 or a pole */
            return (code & 1) ? nfloat_nan(rt, ctx) : nfloat_zero(rt, ctx);
        }

        if (w == NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 3, ctx);

        return _nfloat_set_mpn_err_recip(rt, ys, w + 1, -FLINT_BITS * w, err, recip, code ^ xsgn, ctx);
    }
}

int
nfloat_sin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(res, NULL, NULL, x, ctx);
}

int
nfloat_cos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(NULL, res, NULL, x, ctx);
}

int
nfloat_sin_cos_pi(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res1, ctx) | nfloat_one(res2, ctx);
        return nfloat_nan(res1, ctx) | nfloat_nan(res2, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_trig_pi(res1, res2, NULL, t, ctx);
    }

    return _nfloat_trig_pi(res1, res2, NULL, x, ctx);
}

int
nfloat_tan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(NULL, NULL, res, x, ctx);
}
