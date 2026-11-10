/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"

/* err + c, saturating: a kernel count of UWORD_MAX (no bound, from the
   quotient paths) stays UWORD_MAX */
FLINT_FORCE_INLINE ulong
_err_add(ulong err, ulong c)
{
    return (err > UWORD_MAX - c) ? UWORD_MAX : err + c;
}

/*
    Trigonometric functions of floating-point mantissas (mp_real.h), the
    kernels of the sine, cosine and tangent of nfloat, arf and arb.

    For |x| >= 1, |x| is reduced modulo pi/2 in fixed point by the
    reduction of mp_real_sin_cos_bits (_mp_real_trig_reduce, on the
    mantissa realigned to the limb radix): |x| = a pi/2 + s v, s = +-1,
    v in [0, pi/4] at w fraction limbs within 3 ulps, and

        sin |x| = s sin v, cos v, -s sin v, -cos v   for a = 0, 1, 2, 3,
        cos |x| = cos v, -s sin v, -cos v, s sin v,
        tan |x| = s tan v (a even), -s / tan v (a odd).

    When the needed output is sin v or tan v and v has z leading zero
    bits, the reduction is repeated with z more bits so that the result
    keeps its relative accuracy, as far as the caller's buffers allow
    (wmax).  Arguments below 1 go to the kernels unreduced.

    sin(pi x), cos(pi x), tan(pi x) reduce |x| exactly by the reduction of
    mp_real_sin_cos_pi_bits (_mp_real_trig_pi_reduce): |x| = a/2 + s t
    with t in [0, 1/4], then v = pi t at enough limbs for t's leading
    zero bits; tiny x give pi x (1 - delta) as a product with the
    tabulated pi.
*/

/* the sign bits of sin |x| and cos |x| from the quadrant a and the sign
   s of v */
static int
_signs(int a, int sneg)
{
    int negs, negc;

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

    return negs | (negc << 1);
}

/* x = (d, n) 2^(e - n FLINT_BITS) with e >= 1 as an exact mp_real view
   reading the limbs of buf (n + 1 limbs) */
static void
_view(mp_real_t v, nn_ptr buf, nn_srcptr d, slong n, slong e)
{
    slong q, size;
    unsigned int r;
    nn_ptr p = buf;

    /* |x| = M 2^(e - FLINT_BITS n) = (M 2^r) B^(q + 1 - (n + 1)) with
       e = FLINT_BITS q + r, 0 <= r < FLINT_BITS */
    q = (e >= 0) ? e / FLINT_BITS : -((-e + FLINT_BITS - 1) / FLINT_BITS);
    r = (unsigned int) (e - FLINT_BITS * q);

    if (r == 0)
    {
        flint_mpn_copyi(p, d, n);
        size = n;
        v->exp = q;
    }
    else
    {
        p[n] = mpn_lshift(p, d, n, r);
        size = n + 1;
        v->exp = q + 1;
    }

    while (p[0] == 0)
    {
        p++;
        size--;
    }

    v->d = p;
    v->alloc = size;
    v->size = size;
    v->negative = 0;
    v->err = 0;
}

/* the fraction limbs for z leading zero bits of the reduced argument
   (relative accuracy), with div the margin for the reciprocal to
   normalise with a right shift (MP_REAL_TRIG_FP_NORM_SHIFT bits of left
   shift at most), at most wmax */
static slong
_limbs(slong pb, slong z, int div, slong wmax)
{
    slong w = _mp_real_elem_limbs(pb + z);

    if (div && z > MP_REAL_TRIG_FP_NORM_SHIFT - 1)
        w = FLINT_MAX(w, _mp_real_elem_limbs(pb) + (z - MP_REAL_TRIG_FP_NORM_SHIFT + 1 + FLINT_BITS - 1) / FLINT_BITS);

    return FLINT_MIN(w, wmax);
}

/* sin v and cos v for (v, w) within eps ulps, into the output slots by
   the quadrant (either output NULL) */
static void
_eval(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v, slong w, ulong eps, int a)
{
    nn_ptr os = (a & 1) ? yc : ys;     /* receives sin v */
    nn_ptr oc = (a & 1) ? ys : yc;     /* receives cos v */

    /* (either output NULL: the series paths skip it) */
    _mp_real_sin_cos_kernel(os, oc, err, v, w);
    *err += eps;
}

int
_mp_real_sin_cos_fp(nn_ptr ys, nn_ptr yc, slong * w, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax, int rel)
{
    nn_ptr v, xbuf;
    mp_real_t X;
    slong ww, z, nk = _mp_real_elem_limbs(pb);
    int a, sneg, trunc, r;
    TMP_INIT;

    FLINT_ASSERT(nk <= wmax);
    FLINT_ASSERT(ys != NULL || yc != NULL);
    FLINT_ASSERT(n >= 1 && (d[n - 1] >> (FLINT_BITS - 1)) == 1);
    FLINT_ASSERT(e >= -(WORD_MAX / 2) && e <= MP_REAL_TRIG_MAGLIM(pb));

    if (e <= 0)
    {
        /* |x| < 1 unreduced: z = -e leading zero bits */
        ww = ((rel & 1) && ys != NULL) ? _limbs(pb, -e, 0, wmax) : nk;
        TMP_START;
        v = TMP_ALLOC((ww + 1) * sizeof(ulong));
        trunc = _mp_real_fp_get_fixed(v, ww, ww, d, n, e);
        _eval(ys, yc, err, v, ww, trunc, 0);
        *w = ww;
        TMP_END;
        return 0;
    }

    TMP_START;
    xbuf = TMP_ALLOC((n + 1) * sizeof(ulong));
    _view(X, xbuf, d, n, e);
    v = TMP_ALLOC((wmax + 1) * sizeof(ulong));
    ww = nk;

    for (;;)
    {
        int need_rel;

        r = _mp_real_trig_reduce(v, X, ww);
        a = r & 3;
        sneg = (r & 4) != 0;

        /* the output sin v needs relative accuracy */
        need_rel = ((rel & 1) && ys != NULL && !(a & 1)) || ((rel & 2) && yc != NULL && (a & 1));
        if (!need_rel || ww == wmax)
            break;

        z = _mp_real_elem_lzb(v, ww);
        /* v within 3 ulps: z is accurate unless v is that small */
        if (z >= FLINT_BITS * ww - 4)
            z = FLINT_BITS * ww;

        if (_limbs(pb, z, 0, wmax) <= ww)
            break;
        ww = _limbs(pb, z, 0, wmax);
    }

    _eval(ys, yc, err, v, ww, 3, a);
    *w = ww;
    TMP_END;
    return _signs(a, sneg);
}

int
_mp_real_tan_fp(nn_ptr y, slong * w, ulong * err, int * recip, nn_srcptr d, slong n, slong e, slong pb, slong wmax)
{
    nn_ptr v, xbuf;
    mp_real_t X;
    slong ww, z, nk = _mp_real_elem_limbs(pb);
    int a, sneg, trunc, r;
    TMP_INIT;

    FLINT_ASSERT(nk <= wmax);
    FLINT_ASSERT(n >= 1 && (d[n - 1] >> (FLINT_BITS - 1)) == 1);
    FLINT_ASSERT(e >= -(WORD_MAX / 2) && e <= MP_REAL_TRIG_MAGLIM(pb));

    if (e <= 0)
    {
        /* |x| < 1 unreduced (tan' < 4 there) */
        ww = _limbs(pb, -e, 0, wmax);
        TMP_START;
        v = TMP_ALLOC((ww + 1) * sizeof(ulong));
        trunc = _mp_real_fp_get_fixed(v, ww, ww, d, n, e);
        _mp_real_tan_kernel(y, err, v, ww);
        *err = _err_add(*err, 4 * trunc);
        *w = ww;
        *recip = 0;
        TMP_END;
        return 0;
    }

    TMP_START;
    xbuf = TMP_ALLOC((n + 1) * sizeof(ulong));
    _view(X, xbuf, d, n, e);
    v = TMP_ALLOC((wmax + 1) * sizeof(ulong));
    ww = nk;

    for (;;)
    {
        r = _mp_real_trig_reduce(v, X, ww);
        a = r & 3;
        sneg = (r & 4) != 0;

        if (ww == wmax)
            break;

        z = _mp_real_elem_lzb(v, ww);
        if (z >= FLINT_BITS * ww - 4)
            z = FLINT_BITS * ww;

        if (_limbs(pb, z, a & 1, wmax) <= ww)
            break;
        ww = _limbs(pb, z, a & 1, wmax);
    }

    /* tan v within err ulps, v within 3 (tan' <= 2 on [0, pi/4]) */
    _mp_real_tan_kernel(y, err, v, ww);
    *err = _err_add(*err, 6);
    *w = ww;
    *recip = a & 1;
    TMP_END;

    /* tan |x| = s tan v (a even), -s / tan v (a odd) */
    return (a & 1) ? !sneg : sneg;
}

/* pi |x| for tiny x as a product with pi at pn = max(n, nk) + 2 limbs:
   pi |x| = (y, pn + 1) 2^k within 5 units of 2^k, k = e - FLINT_BITS pn
   (the floor of pi at pn limbs times 4 errs by less than 4 B^-pn, times
   D < B^n; the n low limbs of the product are dropped) */
static void
_pi_times(nn_ptr y, slong * k, ulong * err, nn_srcptr d, slong n, slong e, slong nk)
{
    nn_ptr p, t;
    slong pn = FLINT_MAX(n, nk) + 2;
    TMP_INIT;

    TMP_START;
    p = TMP_ALLOC((pn + 1 + pn + n + 1) * sizeof(ulong));
    t = p + pn + 1;
    p[pn] = mpn_lshift(p, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, pn), pn, 2);
    flint_mpn_mul(t, p, pn + 1, d, n);
    _flint_mpn_copy_small(y, t + n, pn + 1);
    *k = e - FLINT_BITS * pn;
    *err = 5;
    TMP_END;
}

/* whether the tiny path applies: delta < (pi x)^2 < 2^(2e + 4) below a
   unit of the product's frame of pn + 1 limbs */
static int
_tiny(slong e, slong n, slong nk)
{
    slong pn = FLINT_MAX(n, nk) + 2;
    return 2 * e + 4 <= -FLINT_BITS * (pn + 1);
}

/* whether the reduced t = (T, tl) B^(te - tl) in (0, 1/4] is m / 2^j with
   m odd and j <= FLINT_BITS - 8 (a denominator the algebraic methods of
   sin_cos_pi.c take) */
static int
_dyadic_t(ulong * m, slong * j, nn_srcptr T, slong tl, slong te)
{
    unsigned int c;

    if (tl != 1 || te != 0)
        return 0;
    c = flint_ctz(T[0]);
    *j = FLINT_BITS - c;
    *m = T[0] >> c;
    return *j <= FLINT_BITS - 8;
}

int
_mp_real_sin_cos_pi_fp(nn_ptr ys, nn_ptr yc, slong * w, slong * k, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax, int rel)
{
    nn_ptr T, v, xbuf;
    mp_real_t X;
    slong tl, te, z, ww, tj, nk = _mp_real_elem_limbs(pb);
    ulong tm;
    int a, sneg, r;
    TMP_INIT;

    FLINT_ASSERT(nk <= wmax);
    FLINT_ASSERT(ys != NULL || yc != NULL);
    FLINT_ASSERT(n >= 1 && (d[n - 1] >> (FLINT_BITS - 1)) == 1);
    FLINT_ASSERT(e >= -(WORD_MAX / 2));
    *k = 0;

    /* an even integer (all bits of weight >= 2): sin 0, cos 1 */
    if (e > FLINT_BITS * n)
        return MP_REAL_TRIG_EXACT;

    /* tiny x: sin(pi x) = pi x (1 - delta), cos(pi x) = 1 - delta',
       delta < (pi x)^2 < 2^(2e + 4): the product pi x for the sine, with
       delta below one more unit, and 1 for the cosine */
    if (_tiny(e, n, nk))
    {
        if (ys != NULL)
        {
            _pi_times(ys, k, err, d, n, e, nk);
            *err += 1;
        }
        else
            *err = 0;
        if (yc != NULL)
        {
            _flint_mpn_zero_small(yc, nk);
            yc[nk] = 1;
            *err = FLINT_MAX(*err, 1);
        }
        *w = nk;
        return MP_REAL_TRIG_TINY;
    }

    /* |x| = a/2 + s t exactly, t in [0, 1/4] */
    TMP_START;
    xbuf = TMP_ALLOC((n + 1) * sizeof(ulong));
    T = TMP_ALLOC((n + 1) * sizeof(ulong));
    _view(X, xbuf, d, n, e);
    r = _mp_real_trig_pi_reduce(T, &tl, &te, X);
    a = r & 3;
    sneg = (r & 4) != 0;

    if (tl == 0)
    {
        /* a multiple of 1/2: sin, cos in {0, 1, -1} */
        TMP_END;
        return MP_REAL_TRIG_EXACT | a;
    }

    /* 2^-(z + 1) <= t < 2^-z: relative accuracy for t with z leading zero
       bits when the sine of v is wanted */
    z = (slong) flint_clz(T[tl - 1]) - FLINT_BITS * te;
    if (((rel & 1) && ys != NULL && !(a & 1)) || ((rel & 2) && yc != NULL && (a & 1)))
        ww = _limbs(pb, z, 0, wmax);
    else
        ww = nk;

    /* a dyadic t = m / 2^j with a small denominator: the algebraic
       methods of sin_cos_pi.c where they beat the kernel */
    if (_dyadic_t(&tm, &tj, T, tl, te)
        && _mp_real_sin_cos_pi_ui_div_ui_algebraic(UWORD(1) << tj, ww))
    {
        _mp_real_sin_cos_pi_ui_div_ui((a & 1) ? yc : ys, (a & 1) ? ys : yc, err, tm, UWORD(1) << tj, ww);
        *w = ww;
        TMP_END;
        return _signs(a, sneg);
    }

    /* v = pi t within 3 ulps */
    v = TMP_ALLOC((ww + 1) * sizeof(ulong));
    _mp_real_trig_pi_v_fixed(v, T, tl, te, ww);
    v[ww] = 0;

    _eval(ys, yc, err, v, ww, 3, a);
    *w = ww;
    TMP_END;
    return _signs(a, sneg);
}

int
_mp_real_tan_pi_fp(nn_ptr y, slong * w, slong * k, ulong * err, int * recip, nn_srcptr d, slong n, slong e, slong pb, slong wmax)
{
    nn_ptr T, v, xbuf;
    mp_real_t X;
    slong tl, te, z, ww, tj, nk = _mp_real_elem_limbs(pb);
    ulong tm;
    int a, sneg, r;
    TMP_INIT;

    FLINT_ASSERT(nk <= wmax);
    FLINT_ASSERT(n >= 1 && (d[n - 1] >> (FLINT_BITS - 1)) == 1);
    FLINT_ASSERT(e >= -(WORD_MAX / 2));
    *k = 0;
    *recip = 0;

    if (e > FLINT_BITS * n)
        return MP_REAL_TRIG_EXACT;

    /* tiny x: tan(pi x) = pi x (1 + delta) */
    if (_tiny(e, n, nk))
    {
        _pi_times(y, k, err, d, n, e, nk);
        *err += 1;
        *w = nk;
        return MP_REAL_TRIG_TINY;
    }

    TMP_START;
    xbuf = TMP_ALLOC((n + 1) * sizeof(ulong));
    T = TMP_ALLOC((n + 1) * sizeof(ulong));
    _view(X, xbuf, d, n, e);
    r = _mp_real_trig_pi_reduce(T, &tl, &te, X);
    a = r & 3;
    sneg = (r & 4) != 0;

    if (tl == 0)
    {
        /* a multiple of 1/2: 0 (a even) or a pole (a odd) */
        TMP_END;
        return MP_REAL_TRIG_EXACT | a;
    }

    /* t = 1/4: tan(pi/4) = 1, with the sign of the tangent */
    if (tl == 1 && te == 0 && T[0] == (UWORD(1) << (FLINT_BITS - 2)))
    {
        TMP_END;
        return MP_REAL_TRIG_EXACT | MP_REAL_TRIG_ONE | ((a & 1) ? !sneg : sneg);
    }

    z = (slong) flint_clz(T[tl - 1]) - FLINT_BITS * te;
    ww = _limbs(pb, z, a & 1, wmax);

    /* a dyadic t with a small denominator (see above) */
    if (_dyadic_t(&tm, &tj, T, tl, te)
        && _mp_real_tan_pi_ui_div_ui_algebraic(UWORD(1) << tj, ww)
        && (*err = _mp_real_tan_pi_ui_div_ui(y, tm, UWORD(1) << tj, ww)) != UWORD_MAX)
    {
        *w = ww;
        *recip = a & 1;
        TMP_END;
        return (a & 1) ? !sneg : sneg;
    }

    v = TMP_ALLOC((ww + 1) * sizeof(ulong));
    _mp_real_trig_pi_v_fixed(v, T, tl, te, ww);
    v[ww] = 0;

    _mp_real_tan_kernel(y, err, v, ww);
    *err = _err_add(*err, 6);
    *w = ww;
    *recip = a & 1;
    TMP_END;

    return (a & 1) ? !sneg : sneg;
}
