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

/* fraction limbs for a reduced argument with z leading zero bits; with
   div (the cotangent of the reduced argument), enough that the tangent
   normalizes to _nfloat_elem_limbs(ctx, 0) limbs with a left shift of at
   most NFLOAT_ELEM_MAX_NORM_SHIFT bits */
static slong
_nfloat_trig_limbs(gr_ctx_t ctx, slong z, int div)
{
    slong w = _nfloat_elem_limbs(ctx, z);

    if (div && z > NFLOAT_ELEM_MAX_NORM_SHIFT - 1)
        w = FLINT_MAX(w, _nfloat_elem_limbs(ctx, 0) + (z - NFLOAT_ELEM_MAX_NORM_SHIFT + 1 + FLINT_BITS - 1) / FLINT_BITS);

    return w;
}

/* the trigonometric functions of v (sin and cos for which = 0, tan for
   which = 1) with |x| = a pi/2 + s v (or pi |x| = ...), assembled into
   rs, rc (either may be NULL) resp. rt; v at w limbs within eps ulps */
static int
_nfloat_trig_eval(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nn_srcptr v, slong w,
    ulong eps, int a, int s, int xsgn, gr_ctx_t ctx)
{
    ulong ys[NFLOAT_ELEM_MAX_LIMBS + 2], yc[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    int status = GR_SUCCESS;

    if (rt == NULL)
    {
        int sneg = (s < 0), negs, negc;
        nfloat_ptr os = (a & 1) ? rc : rs;
        nfloat_ptr oc = (a & 1) ? rs : rc;

        _mp_real_sin_cos_kernel(ys, yc, &err, v, w);
        err += eps;

        /* sin |x| = s sin v, cos v, -s sin v, -cos v;
           cos |x| = cos v, -s sin v, -cos v, s sin v   (a = 0, 1, 2, 3) */
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

        /* os gets sin v, oc gets cos v; the sine of x is odd */
        if (os != NULL)
            status |= _nfloat_set_mpn_err(os, ys, w + 1, -FLINT_BITS * w, err, err,
                ((a & 1) ? negc : negs) ^ ((os == rs) ? xsgn : 0), ctx);
        if (oc != NULL)
            status |= _nfloat_set_mpn_err(oc, yc, w + 1, -FLINT_BITS * w, err, err,
                ((a & 1) ? negs : negc) ^ ((oc == rs) ? xsgn : 0), ctx);
        return status;
    }
    else
    {
        int neg;

        _mp_real_tan_kernel(ys, &err, v, w);
        err += 2 * eps;     /* tan' <= 2 on [0, pi/4] */

        if (!(a & 1))
        {
            /* s tan v */
            neg = (s < 0) ^ xsgn;
            return _nfloat_set_mpn_err(rt, ys, w + 1, -FLINT_BITS * w, err, err, neg, ctx);
        }
        else
        {
            /* -s / tan v: tan v normalized to wn limbs, then one division */
            ulong t[NFLOAT_ELEM_MAX_LIMBS + 1], one[NFLOAT_ELEM_MAX_LIMBS + 1];
            ulong q[NFLOAT_ELEM_MAX_LIMBS + 2];
            slong te, wn = _nfloat_elem_limbs(ctx, 0);
            ulong et, eq;

            neg = (s > 0) ^ xsgn;

            if (flint_mpn_zero_p(ys, w + 1))
                return GR_UNABLE;

            et = _nfloat_elem_normalize(t, &te, ys, w + 1, -FLINT_BITS * w, err, wn);
            if (et == UWORD_MAX)
                return GR_UNABLE;
            _nfloat_zero_limbs(one, wn);
            one[wn - 1] = UWORD(1) << (FLINT_BITS - 1);
            eq = _nfloat_elem_div(q, one, 0, t, et, wn);

            /* 1 / (T 2^te) = Q 2^(1 - 2 FLINT_BITS wn - te) */
            return _nfloat_set_mpn_err(rt, q, wn + 1, 1 - 2 * FLINT_BITS * wn - te, eq, eq, neg, ctx);
        }
    }
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
#define NFLOAT_TRIG_FAST_SIGNS(negs_out, negc_out, a, sneg, xsgn) \
    do { \
        int __negs, __negc; \
        if ((a) & 1) { __negs = ((a) == 3); __negc = (sneg) ^ ((a) == 1); } \
        else { __negs = (sneg) ^ ((a) == 2); __negc = ((a) == 2); } \
        /* sin |x| from sin v (a even) or cos v (a odd), and conversely */ \
        negs_out = ((a) & 1) ? __negc : __negs; \
        negc_out = ((a) & 1) ? __negs : __negc; \
    } while (0)

/* reduction of x: v = (v1, v0) at two fraction limbs (ev = 2 for |x| >= 1) */
FLINT_FORCE_INLINE int
_nfloat_trig_reduce_fast_1(int * ra, int * rsneg, ulong * v1, ulong * v0, ulong d, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
    ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
    ulong xi, x1, q, q2, h, l, h2, l2, p3, p2, p1, p0, t3, t2, t1, t0, u3, u2, u1, u0;
    int sneg = 0;

    /* |x| = (xi, x1) with the binary point after xi; 0 <= e < FLINT_BITS */
    xi = (e == 0) ? 0 : (d >> (FLINT_BITS - e));
    x1 = (e == 0) ? d : (d << e);

    /* q = floor(|x| 2/pi), or one less (the one-limb K errs by
       |x| 2^-64 < 1/2) */
    umul_ppmm(h, l, x1, K);
    umul_ppmm(h2, l2, xi, K);
    add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
    (void) l;

    /* t = |x| - 2q P >= 0 */
    q2 = 2 * q;
    umul_ppmm(p1, p0, q2, P[0]);
    umul_ppmm(h, l, q2, P[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q2, P[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(t3, t2, t1, t0, xi, x1, UWORD(0), UWORD(0), p3, p2, p1, p0);

    /* u = 2P (pi/2) */
    u3 = P[2] >> (FLINT_BITS - 1);
    u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
    u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    u0 = P[0] << 1;

    if (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2 && (t1 > u1 || (t1 == u1 && t0 >= u0))))))
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, t3, t2, t1, t0, u3, u2, u1, u0);
        q++;
    }

    /* fold t > pi/4 */
    if (t3 != 0 || t2 > P[2] || (t2 == P[2] && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, u3, u2, u1, u0, t3, t2, t1, t0);
        q++;
        sneg = 1;
    }

    if (t2 < (UWORD(1) << 32))
        return NFLOAT_TRIG_FAST_DECLINE;

    *ra = q & 3;
    *rsneg = sneg;
    *v1 = t2;
    *v0 = t1;
    return 0;
}

FLINT_FORCE_INLINE int
_nfloat_trig_reduce_fast_2(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong d1, ulong d0, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
    ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
    ulong xi, x2, x1, q, q2, h, l, h2, l2, p4, p3, p2, p1, p0;
    ulong t4, t3, t2, t1, t0, u4, u3, u2, u1, u0;
    int sneg = 0;

    /* |x| = (xi, x2, x1) with the binary point after xi; 0 <= e < FLINT_BITS */
    if (e == 0)
    {
        xi = 0; x2 = d1; x1 = d0;
    }
    else
    {
        xi = d1 >> (FLINT_BITS - e);
        x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        x1 = d0 << e;
    }

    umul_ppmm(h, l, x2, K);
    umul_ppmm(h2, l2, xi, K);
    add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
    (void) l;

    q2 = 2 * q;
    umul_ppmm(p1, p0, q2, P[0]);
    umul_ppmm(h, l, q2, P[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q2, P[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, q2, P[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, xi, x2, x1, UWORD(0), UWORD(0), p4, p3, p2, p1, p0);

    u4 = P[3] >> (FLINT_BITS - 1);
    u3 = (P[3] << 1) | (P[2] >> (FLINT_BITS - 1));
    u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
    u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    u0 = P[0] << 1;

    if (t4 > u4 || (t4 == u4 && (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2
            && (t1 > u1 || (t1 == u1 && t0 >= u0))))))))
    {
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, t4, t3, t2, t1, t0, u4, u3, u2, u1, u0);
        q++;
    }

    if (t4 != 0 || t3 > P[3] || (t3 == P[3] && (t2 > P[2] || (t2 == P[2]
            && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))))
    {
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, u4, u3, u2, u1, u0, t4, t3, t2, t1, t0);
        q++;
        sneg = 1;
    }

    if (t3 < (UWORD(1) << 24))
        return NFLOAT_TRIG_FAST_DECLINE;

    *ra = q & 3;
    *rsneg = sneg;
    *v2 = t3;
    *v1 = t2;
    *v0 = t1;
    return 0;
}

/* v = pi t = 3t + (pi - 3) t for t = (t2, t1, t0) <= 1/4 */
FLINT_FORCE_INLINE void
_nfloat_trig_pi_times(ulong * v2, ulong * v1, ulong * v0, ulong t2, ulong t1, ulong t0)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
    ulong c2, c1, c0, s2, s1, s0;

    /* pi - 3 at three fraction limbs, floored: the fraction of 4P */
    c2 = (P[3] << 2) | (P[2] >> (FLINT_BITS - 2));
    c1 = (P[2] << 2) | (P[1] >> (FLINT_BITS - 2));
    c0 = (P[1] << 2) | (P[0] >> (FLINT_BITS - 2));

    _mp_real_small_mulhi_3x3(&s2, &s1, &s0, c2, c1, c0, t2, t1, t0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, t2, t1, t0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, (t2 << 1) | (t1 >> (FLINT_BITS - 1)),
        (t1 << 1) | (t0 >> (FLINT_BITS - 1)), t0 << 1);
    *v2 = s2; *v1 = s1; *v0 = s0;
}

/* 2|x| = m + (f2, f1, f0) to the reduction pi |x| = a pi/2 + s v: t at
   three limbs (exact if f0's lowest bit is zero) */
FLINT_FORCE_INLINE int
_nfloat_trig_pi_fold(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong m, ulong f2, ulong f1, ulong f0)
{
    if (f2 == 0 && f1 == 0 && f0 == 0)
        return NFLOAT_TRIG_FAST_DECLINE;

    if (f2 < (UWORD(1) << (FLINT_BITS - 1)) || (f2 == (UWORD(1) << (FLINT_BITS - 1)) && f1 == 0 && f0 == 0))
    {
        *ra = m & 3;
        *rsneg = 0;
    }
    else
    {
        *ra = (m + 1) & 3;
        *rsneg = 1;
        sub_dddmmmsss(f2, f1, f0, UWORD(0), UWORD(0), UWORD(0), f2, f1, f0);
    }

    /* t = f / 2 */
    _nfloat_trig_pi_times(v2, v1, v0, f2 >> 1, (f2 << (FLINT_BITS - 1)) | (f1 >> 1),
        (f1 << (FLINT_BITS - 1)) | (f0 >> 1));
    return 0;
}

/* reduction of pi x, one limb, -31 <= e < FLINT_BITS */
FLINT_FORCE_INLINE int
_nfloat_trig_pi_reduce_fast_1(int * ra, int * rsneg, ulong * v1, ulong * v0, ulong d, slong e)
{
    ulong m, f1, f0, v2, w1, w0;

    if (e >= 1)
    {
        if (e == FLINT_BITS - 1)
            return NFLOAT_TRIG_FAST_DECLINE;
        m = d >> (FLINT_BITS - 1 - e);
        f1 = d << (e + 1);
        f0 = 0;
    }
    else if (e == 0)
    {
        m = 1;
        f1 = d << 1;
        f0 = 0;
    }
    else
    {
        unsigned int s = -e - 1;
        m = 0;
        f1 = d >> s;
        f0 = (s == 0) ? 0 : (d << (FLINT_BITS - s));
    }

    if (_nfloat_trig_pi_fold(ra, rsneg, &v2, &w1, &w0, m, f1, f0, 0) != 0)
        return NFLOAT_TRIG_FAST_DECLINE;

    if (v2 < (UWORD(1) << 32))
        return NFLOAT_TRIG_FAST_DECLINE;

    *v1 = v2;
    *v0 = w1;
    return 0;
}

/* reduction of pi x, two limbs, -39 <= e < FLINT_BITS */
FLINT_FORCE_INLINE int
_nfloat_trig_pi_reduce_fast_2(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong d1, ulong d0, slong e)
{
    ulong m, f2, f1, f0;

    if (e >= 1)
    {
        if (e == FLINT_BITS - 1)
        {
            m = d1;
            f2 = d0;
            f1 = 0;
        }
        else
        {
            m = d1 >> (FLINT_BITS - 1 - e);
            f2 = (d1 << (e + 1)) | (d0 >> (FLINT_BITS - 1 - e));
            f1 = d0 << (e + 1);
        }
        f0 = 0;
    }
    else if (e == 0)
    {
        m = 1;
        f2 = (d1 << 1) | (d0 >> (FLINT_BITS - 1));
        f1 = d0 << 1;
        f0 = 0;
    }
    else
    {
        unsigned int s = -e - 1;
        m = 0;
        if (s == 0)
        {
            f2 = d1; f1 = d0; f0 = 0;
        }
        else
        {
            f2 = d1 >> s;
            f1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
            f0 = d0 << (FLINT_BITS - s);
        }
    }

    if (_nfloat_trig_pi_fold(ra, rsneg, v2, v1, v0, m, f2, f1, f0) != 0)
        return NFLOAT_TRIG_FAST_DECLINE;

    if (*v2 < (UWORD(1) << 24))
        return NFLOAT_TRIG_FAST_DECLINE;

    return 0;
}

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

            NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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

            NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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

    NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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

            NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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

            NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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

    NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
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
            if (e < -31 || _nfloat_trig_pi_reduce_fast_1(&a, &sneg, &v1, &v0, d, e) != 0)
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
            if (_nfloat_trig_reduce_fast_1(&a, &sneg, &v1, &v0, d, e) != 0)
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
            if (e < -39 || _nfloat_trig_pi_reduce_fast_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
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
            if (_nfloat_trig_reduce_fast_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
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
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 1];
    ulong xbuf[NFLOAT_MAX_LIMBS + 1];
    mp_real_t X;
    slong n, e, w, z;
    int xsgn, a, s;
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

    if (e <= 0)
    {
        /* |x| < 1 unreduced: v = |x| with z = -e leading zero bits */
        int trunc;
        z = -e;
        w = _nfloat_elem_limbs(ctx, (rs != NULL || rt != NULL) ? z : 0);
        trunc = _nfloat_get_fixed(v, w, w, NFLOAT_D(x), n, e);
        return _nfloat_trig_eval(rs, rc, rt, v, w, trunc, 0, 1, xsgn, ctx);
    }

    if (e > NFLOAT_TRIG_MAX_EXP)
        return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);

    w = _nfloat_elem_limbs(ctx, 0);
    _nfloat_mp_real_view(X, xbuf, x, ctx);

    for (;;)
    {
        int need_rel, r;

        /* the reduction of mp_real_sin_cos_bits */
        r = _mp_real_trig_reduce(v, X, w);
        a = r & 3;
        s = (r & 4) ? -1 : 1;

        /* the output sin v (or tan v) needs relative accuracy */
        need_rel = (rt != NULL) || (rs != NULL && !(a & 1)) || (rc != NULL && (a & 1));

        if (!need_rel)
            break;

        z = _mp_real_elem_lzb(v, w);

        /* v within 3 ulps: z is accurate unless v is that small */
        if (z >= FLINT_BITS * w - 4)
            z = FLINT_BITS * w;

        if (_nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1)) <= w)
            break;

        w = _nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1));

        /* more cancellation than the temporaries allow */
        if (w > NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);
    }

    return _nfloat_trig_eval(rs, rc, rt, v, w, 3, a, s, xsgn, ctx);
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
    ulong T[NFLOAT_MAX_LIMBS + 1], xbuf[NFLOAT_MAX_LIMBS + 1];
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 2];
    mp_real_t X;
    slong n, e, w, z, tl, te;
    int xsgn, a, s, r;
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

    /* |x| = a/2 + s t exactly, t in [0, 1/4] (the reduction of
       mp_real_sin_cos_pi_bits) */
    _nfloat_mp_real_view(X, xbuf, x, ctx);
    r = _mp_real_trig_pi_reduce(T, &tl, &te, X);
    a = r & 3;
    s = (r & 4) ? -1 : 1;

    if (tl == 0)
    {
        /* x a multiple of 1/2: sin, cos in {0, 1, -1}, tan in {0, pole} */
        if (rt != NULL)
            return (a & 1) ? nfloat_nan(rt, ctx) : nfloat_zero(rt, ctx);

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

    /* 2^-(z + 1) <= t < 2^-z: relative accuracy for t with z leading
       zero bits */
    z = (slong) flint_clz(T[tl - 1]) - FLINT_BITS * te;
    w = _nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1));
    if (w > NFLOAT_ELEM_MAX_LIMBS)
        return _nfloat_trig_slow(rs, rc, rt, x, 3, ctx);

    /* v = pi t within 3 ulps */
    _mp_real_trig_pi_v_fixed(v, T, tl, te, w);
    v[w] = 0;

    return _nfloat_trig_eval(rs, rc, rt, v, w, 3, a, s, xsgn, ctx);
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
