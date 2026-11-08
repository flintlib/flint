/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The real factors of the complex elementary functions at one and two
   limbs (the kernels of acf/elem.h compiled once, for acf, nfloat_complex,
   acb and arb). */

#include "elem.h"

#if FLINT_BITS == 64

/* |sin x| = S 2^(ks - 2 FLINT_BITS) and |cos x| = C 2^-(2 FLINT_BITS) for
   the frames S = (s2, s1, s0) and C = (c2, c1, c0) (values below 1)
   within *err units, with the signs *negs of sin x and *negc of cos x,
   for x = (-1)^xsgn d 2^(e - FLINT_BITS), or of pi x (pi = 1).  Tiny
   arguments (e <= _ACF_SMALL_TINY_EXP_1) give sin x = x and cos x = 1
   within two units (the frame of x, or of pi |x| with its lowest limb
   dropped, has ks = e).  Returns 0, or ARF_ELEM_UNCERTAIN when the
   reduction declines. */
int
_acf_small_trig_1(ulong * s2, ulong * s1, ulong * s0, slong * ks, int * negs,
    ulong * c2, ulong * c1, ulong * c0, int * negc, ulong * err, ulong d, slong e, int xsgn, int pi)
{
    ulong v1, v0, t1, t0, g1, g0, ev;
    int a = 0, sneg = 0, ns, nc;

    if (e >= FLINT_BITS)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_1)
    {
        if (pi)
        {
            /* pi (d / B) with the lowest limb dropped */
            ulong y0;
            _acf_small_pi_times_1(s2, s1, s0, &y0, d, 0);
            *err = 4;
        }
        else
        {
            *s2 = 0; *s1 = d; *s0 = 0;
            *err = 2;
        }
        *ks = e;
        *c2 = 1; *c1 = 0; *c0 = 0;
        *negs = xsgn;
        *negc = 0;
        return 0;
    }

    if (pi)
    {
        if (e < -31)
        {
            /* v = pi |x| < 1 directly (the lowest limb of the product
               dropped) */
            ulong xi, x0;
            _acf_small_pi_times_1(&xi, &v1, &v0, &x0, d, e);
            ev = 2;
        }
        else if (_arf_small_trig_pi_reduce_1(&a, &sneg, &v1, &v0, d, e) != 0)
            return ARF_ELEM_UNCERTAIN;
        else
            ev = 2;
    }
    else if (e <= 0)
    {
        /* v = |x| at two fraction limbs (truncated below 2^-(2 FLINT_BITS)) */
        ulong t;
        _arf_small_rshift_3(&t, &v1, &v0, UWORD(0), d, UWORD(0), -e);
        ev = (e < -FLINT_BITS);
    }
    else
    {
        if (_arf_small_trig_reduce_1(&a, &sneg, &v1, &v0, d, e) != 0)
            return ARF_ELEM_UNCERTAIN;
        ev = 2;
    }

    _mp_real_small_sin_cos_g_2(&t1, &t0, &g1, &g0, v1, v0);
    *err = MP_REAL_SMALL_SIN_COS_ERR_2 + ev;
    /* (the macro gives the signs carried by the frames of sin v and
       cos v, which feed sin |x| and cos |x| for even a, cos |x| and
       sin |x| for odd a) */
    ARF_SMALL_TRIG_SIGNS(ns, nc, a, sneg, xsgn);
    *negs = ((a & 1) ? nc : ns) ^ xsgn;
    *negc = (a & 1) ? ns : nc;
    *ks = 0;

    if (a & 1)
    {
        sub_dddmmmsss(*s2, *s1, *s0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
        *c2 = 0; *c1 = t1; *c0 = t0;
    }
    else
    {
        *s2 = 0; *s1 = t1; *s0 = t0;
        sub_dddmmmsss(*c2, *c1, *c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
    }
    return 0;
}

/* the same at two limbs for x = (-1)^xsgn (d1, d0) 2^(e - 2 FLINT_BITS):
   frames with three fraction limbs */
int
_acf_small_trig_2(ulong * s3, ulong * s2, ulong * s1, ulong * s0, slong * ks, int * negs,
    ulong * c3, ulong * c2, ulong * c1, ulong * c0, int * negc, ulong * err,
    ulong d1, ulong d0, slong e, int xsgn, int pi)
{
    ulong v2, v1, v0, t2, t1, t0, g2, g1, g0, ev;
    int a = 0, sneg = 0, ns, nc;

    if (e >= FLINT_BITS)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_2)
    {
        if (pi)
        {
            ulong y0;
            _acf_small_pi_times_2(s3, s2, s1, s0, &y0, d1, d0, 0);
            *err = 4;
        }
        else
        {
            *s3 = 0; *s2 = d1; *s1 = d0; *s0 = 0;
            *err = 2;
        }
        *ks = e;
        *c3 = 1; *c2 = 0; *c1 = 0; *c0 = 0;
        *negs = xsgn;
        *negc = 0;
        return 0;
    }

    if (pi)
    {
        if (e < -39)
        {
            ulong xi, x0;
            _acf_small_pi_times_2(&xi, &v2, &v1, &v0, &x0, d1, d0, e);
            ev = 2;
        }
        else if (_arf_small_trig_pi_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
            return ARF_ELEM_UNCERTAIN;
        else
            ev = 8;
    }
    else if (e <= 0)
    {
        _arf_small_rshift_3(&v2, &v1, &v0, d1, d0, UWORD(0), -e);
        ev = (e < -FLINT_BITS);
    }
    else
    {
        if (_arf_small_trig_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
            return ARF_ELEM_UNCERTAIN;
        ev = 2;
    }

    _mp_real_small_sin_cos_g_3(&t2, &t1, &t0, &g2, &g1, &g0, v2, v1, v0);
    *err = MP_REAL_SMALL_SIN_COS_ERR_3 + ev;
    /* (the macro gives the signs carried by the frames of sin v and
       cos v, which feed sin |x| and cos |x| for even a, cos |x| and
       sin |x| for odd a) */
    ARF_SMALL_TRIG_SIGNS(ns, nc, a, sneg, xsgn);
    *negs = ((a & 1) ? nc : ns) ^ xsgn;
    *negc = (a & 1) ? ns : nc;
    *ks = 0;

    if (a & 1)
    {
        sub_ddddmmmmssss(*s3, *s2, *s1, *s0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
        *c3 = 0; *c2 = t2; *c1 = t1; *c0 = t0;
    }
    else
    {
        *s3 = 0; *s2 = t2; *s1 = t1; *s0 = t0;
        sub_ddddmmmmssss(*c3, *c2, *c1, *c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
    }
    return 0;
}

/* |sinh x| = S 2^(ks - 2 FLINT_BITS) and cosh x = C 2^(kc - 2 FLINT_BITS)
   (values below _ACF_SMALL_HYP_BOUND) within *err units for |x| = d
   2^(e - FLINT_BITS), or of pi |x| (pi = 1); tiny arguments as for the
   trigonometric functions; e <= MP_REAL_EXP_FP_MAX_EXP - 2 */
int
_acf_small_hyp_1(ulong * s2, ulong * s1, ulong * s0, slong * ks,
    ulong * c2, ulong * c1, ulong * c0, slong * kc, ulong * err, ulong d, slong e, int pi)
{
    ulong y2, y1, y0, z2, z1, z0, q;

    if (e > MP_REAL_EXP_FP_MAX_EXP - 2)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_1)
    {
        if (pi)
        {
            /* pi (d / B) with the lowest limb dropped */
            ulong y0;
            _acf_small_pi_times_1(s2, s1, s0, &y0, d, 0);
            *err = 4;
        }
        else
        {
            *s2 = 0; *s1 = d; *s0 = 0;
            *err = 2;
        }
        *ks = e;
        *c2 = 1; *c1 = 0; *c0 = 0;
        *kc = 0;
        return 0;
    }

    if (pi)
    {
        ulong xi, x2, x1, x0;
        _acf_small_pi_times_1(&xi, &x2, &x1, &x0, d, e);
        q = _arf_small_exp_pair_fixed_1(&y2, &y1, &y0, &z2, &z1, &z0, xi, x2, x1, x0);
        *err = ARF_SMALL_HYP_ERR_2 + 2;
    }
    else
    {
        q = _arf_small_exp_pair_1(&y2, &y1, &y0, &z2, &z1, &z0, d, e);
        *err = ARF_SMALL_HYP_ERR_2;
    }

    *ks = *kc = (slong) q - 1;
    sub_dddmmmsss(*s2, *s1, *s0, y2, y1, y0, z2, z1, z0);
    add_sssaaaaaa(*c2, *c1, *c0, y2, y1, y0, z2, z1, z0);
    return 0;
}

int
_acf_small_hyp_2(ulong * s3, ulong * s2, ulong * s1, ulong * s0, slong * ks,
    ulong * c3, ulong * c2, ulong * c1, ulong * c0, slong * kc, ulong * err, ulong d1, ulong d0, slong e, int pi)
{
    ulong y3, y2, y1, y0, z3, z2, z1, z0, q;

    if (e > MP_REAL_EXP_FP_MAX_EXP - 2)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_2)
    {
        if (pi)
        {
            ulong y0;
            _acf_small_pi_times_2(s3, s2, s1, s0, &y0, d1, d0, 0);
            *err = 4;
        }
        else
        {
            *s3 = 0; *s2 = d1; *s1 = d0; *s0 = 0;
            *err = 2;
        }
        *ks = e;
        *c3 = 1; *c2 = 0; *c1 = 0; *c0 = 0;
        *kc = 0;
        return 0;
    }

    if (pi)
    {
        ulong xi, x3, x2, x1, x0;
        _acf_small_pi_times_2(&xi, &x3, &x2, &x1, &x0, d1, d0, e);
        q = _arf_small_exp_pair_fixed_2(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, xi, x3, x2, x1, x0);
        *err = ARF_SMALL_HYP_ERR_3 + 2;
    }
    else
    {
        q = _arf_small_exp_pair_2(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, d1, d0, e);
        *err = ARF_SMALL_HYP_ERR_3;
    }

    *ks = *kc = (slong) q - 1;
    sub_ddddmmmmssss(*s3, *s2, *s1, *s0, y3, y2, y1, y0, z3, z2, z1, z0);
    add_ssssaaaaaaaa(*c3, *c2, *c1, *c0, y3, y2, y1, y0, z3, z2, z1, z0);
    return 0;
}

/* E = exp(x) (frame value below 2) for x = (-1)^sgnbit d 2^(e - FLINT_BITS),
   or exp(-pi x) (pi = 1), as (e2, e1, e0) 2^(k - 2 FLINT_BITS) within
   *err units; e <= MP_REAL_EXP_FP_MAX_EXP - 2.  Tiny arguments by the
   linear term (the quadratic one below a unit). */
int
_acf_small_exp_factor_1(ulong * e2, ulong * e1, ulong * e0, slong * k, ulong * err, ulong d, slong e, int sgnbit, int pi)
{
    ulong xi, x2, x1, x0;

    if (e > MP_REAL_EXP_FP_MAX_EXP - 2)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_1)
    {
        /* |x| or pi |x| in units of 2^-(2 FLINT_BITS): the frame of
           |x| / 2^e, shifted down by -e bits */
        if (pi)
        {
            ulong x3;
            _acf_small_pi_times_1(&xi, &x2, &x1, &x0, d, 0);
            _arf_small_rshift_4(&x3, &x2, &x1, &x0, xi, x2, x1, x0, FLINT_BITS - e);
            sgnbit ^= 1;
        }
        else
        {
            x2 = 0; x1 = 0;
            x0 = (e <= -2 * FLINT_BITS) ? 0 : (d >> (-FLINT_BITS - e));
        }
        if (sgnbit)
            sub_dddmmmsss(*e2, *e1, *e0, UWORD(1), UWORD(0), UWORD(0), x2, x1, x0);
        else
            add_sssaaaaaa(*e2, *e1, *e0, UWORD(1), UWORD(0), UWORD(0), x2, x1, x0);
        *k = 0;
        *err = 4;
        return 0;
    }

    if (pi)
    {
        _acf_small_pi_times_1(&xi, &x2, &x1, &x0, d, e);
        *k = _arf_small_exp_fixed_1(e2, e1, e0, xi, x2, x1, x0, !sgnbit);
        *err = ARF_SMALL_EXP_ERR_2 + 1;
    }
    else
    {
        *k = _arf_small_exp_1(e2, e1, e0, d, e, sgnbit);
        *err = ARF_SMALL_EXP_ERR_2;
    }
    return 0;
}

int
_acf_small_exp_factor_2(ulong * e3, ulong * e2, ulong * e1, ulong * e0, slong * k, ulong * err,
    ulong d1, ulong d0, slong e, int sgnbit, int pi)
{
    ulong xi, x3, x2, x1, x0;

    if (e > MP_REAL_EXP_FP_MAX_EXP - 2)
        return ARF_ELEM_UNCERTAIN;

    if (e <= _ACF_SMALL_TINY_EXP_2)
    {
        if (pi)
        {
            /* (the lowest limb of the product dropped first) */
            _acf_small_pi_times_2(&xi, &x3, &x2, &x1, &x0, d1, d0, 0);
            _arf_small_rshift_4(&x3, &x2, &x1, &x0, xi, x3, x2, x1, -e);
            sgnbit ^= 1;
        }
        else
        {
            /* (d1, d0) 2^(e - 2 FLINT_BITS) in units of 2^-(3 FLINT_BITS) */
            x3 = 0;
            _arf_small_rshift_3(&x2, &x1, &x0, UWORD(0), d1, d0, -FLINT_BITS - e);
        }
        if (sgnbit)
            sub_ddddmmmmssss(*e3, *e2, *e1, *e0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), x3, x2, x1, x0);
        else
            add_ssssaaaaaaaa(*e3, *e2, *e1, *e0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), x3, x2, x1, x0);
        *k = 0;
        *err = 4;
        return 0;
    }

    if (pi)
    {
        _acf_small_pi_times_2(&xi, &x3, &x2, &x1, &x0, d1, d0, e);
        *k = _arf_small_exp_fixed_2(e3, e2, e1, e0, xi, x3, x2, x1, x0, !sgnbit);
        *err = ARF_SMALL_EXP_ERR_3 + 1;
    }
    else
    {
        *k = _arf_small_exp_2(e3, e2, e1, e0, d1, d0, e, sgnbit);
        *err = ARF_SMALL_EXP_ERR_3;
    }
    return 0;
}

/* log(a^2 + b^2) / 2 for a = da 2^(ea - FLINT_BITS), b = db 2^(eb -
   FLINT_BITS) (nonzero mantissas) as the frame re: the sum of squares
   exactly (the smaller square shifted into place, its dropped bits
   below a unit), then the logarithm kernel of arf at two limbs, or
   near 1 the series of log(1 + D) on the exact D = a^2 + b^2 - 1 */
int
_acf_small_log_abs_1(_acf_frame_struct * re, ulong da, slong ea, ulong db, slong eb)
{
    ulong h1, l1, h2, l2, s2, s1, s0, m1, m0, r2, r1, r0, derr;
    slong E, sh, X;
    unsigned int c;

    /* S = (s2, s1, s0) 2^(2E - 2 FLINT_BITS), E the larger exponent */
    umul_ppmm(h1, l1, da, da);
    umul_ppmm(h2, l2, db, db);
    if (ea >= eb)
    {
        E = ea;
        sh = 2 * (ea - eb);
        _arf_small_rshift_3(&s2, &h2, &l2, UWORD(0), h2, l2, sh);
    }
    else
    {
        E = eb;
        sh = 2 * (eb - ea);
        _arf_small_rshift_3(&s2, &h1, &l1, UWORD(0), h1, l1, sh);
    }
    add_sssaaaaaa(s2, s1, s0, UWORD(0), h1, l1, UWORD(0), h2, l2);
    derr = (sh > 0);

    /* the normalised two-limb window M 2^X of S 2^(2E - 2 FLINT_BITS)
       (S >= 2^126: s2 <= 1, and the window starts in s2 or s1) */
    if (s2 != 0)
    {
        m1 = (s2 << (FLINT_BITS - 1)) | (s1 >> 1);
        m0 = (s1 << (FLINT_BITS - 1)) | (s0 >> 1);
        X = 2 * E + 1;
        c = FLINT_BITS - 1;
    }
    else
    {
        c = flint_clz(s1);
        m1 = (c == 0) ? s1 : ((s1 << c) | (s0 >> (FLINT_BITS - c)));
        m0 = (c == 0) ? s0 : (s0 << c);
        X = 2 * E - (slong) c;
    }

    if (ARF_SMALL_LOG_OK_1(m1, X))
    {
        /* (the window's and the sum's dropped bits: 2^-126 relative,
           below 4 units of the logarithm's frame) */
        re->sgnbit = _arf_small_log_2(&r2, &r1, &r0, m1, m0, X);
        re->y[2] = r2; re->y[1] = r1; re->y[0] = r0;
        re->k = -1;
        re->err = ARF_SMALL_LOG_ERR_2 + 6;
        return 0;
    }
    else
    {
        /* D = S 2^(2E - 2 FLINT_BITS) - 1 exactly: the unit lies at bit
           2 FLINT_BITS - 2E of S (X = 0 or 1) */
        ulong d1, d0;
        slong p = 2 * FLINT_BITS - 2 * E, ed;
        int dneg;

        if (s2 != 0)
        {
            FLINT_ASSERT(p == 2 * FLINT_BITS);
            d1 = s1; d0 = s0;
            dneg = 0;
        }
        else if (p == 2 * FLINT_BITS)
        {
            /* S < 1: D = -(1 - S) */
            sub_ddmmss(d1, d0, UWORD(0), UWORD(0), s1, s0);
            dneg = 1;
        }
        else
        {
            ulong u1;
            FLINT_ASSERT(p >= FLINT_BITS && p < 2 * FLINT_BITS);
            u1 = UWORD(1) << (p - FLINT_BITS);
            if (s1 >= u1)
            {
                sub_ddmmss(d1, d0, s1, s0, u1, UWORD(0));
                dneg = 0;
            }
            else
            {
                sub_ddmmss(d1, d0, u1, UWORD(0), s1, s0);
                dneg = 1;
            }
        }

        if (d1 == 0 && d0 == 0)
        {
            /* |z| = 1 exactly: log = 0 (the frame of a zero is declined) */
            return ARF_ELEM_UNCERTAIN;
        }

        /* normalise D = (d1, d0) 2^(ed - 2 FLINT_BITS) */
        if (d1 != 0)
        {
            c = flint_clz(d1);
            if (c != 0)
            {
                d1 = (d1 << c) | (d0 >> (FLINT_BITS - c));
                d0 = d0 << c;
            }
            ed = 2 * E - (slong) c;
        }
        else
        {
            c = flint_clz(d0);
            d1 = d0 << c;
            d0 = 0;
            ed = 2 * E - FLINT_BITS - (slong) c;
        }

        if (ed > -31)
            return ARF_ELEM_UNCERTAIN;

        /* the sum's dropped bits, a unit of S's frame, are 2^(2E - ed)
           units of D's */
        if (derr)
        {
            if (2 * E - ed > 56)
                return ARF_ELEM_UNCERTAIN;
            derr = UWORD(1) << (2 * E - ed);
        }

        _acf_small_log1p_series_2(&r2, &r1, &r0, &re->err, d1, d0, ed, dneg);
        re->y[2] = r2; re->y[1] = r1; re->y[0] = r0;
        re->k = ed - 1;
        re->sgnbit = dneg;
        re->err += 2 * derr;
        return 0;
    }
}

int
_acf_small_log_abs_2(_acf_frame_struct * re, ulong da1, ulong da0, slong ea, ulong db1, ulong db0, slong eb)
{
    ulong A[4], Bq[4], s[5], m2, m1, m0, r3, r2, r1, r0, derr, lerr = 0;
    slong E, sh, X;
    unsigned int c;

    /* S = (s4, s3, s2, s1, s0) 2^(2E - 4 FLINT_BITS) */
    {
        ulong t[2];
        t[0] = da0; t[1] = da1;
        flint_mpn_sqr(A, t, 2);
        t[0] = db0; t[1] = db1;
        flint_mpn_sqr(Bq, t, 2);
    }
    if (ea >= eb)
    {
        E = ea;
        sh = 2 * (ea - eb);
        s[4] = 0;
        s[3] = A[3]; s[2] = A[2]; s[1] = A[1]; s[0] = A[0];
        if (sh < 4 * FLINT_BITS)
        {
            ulong t[4];
            slong sl = sh / FLINT_BITS, sb = sh % FLINT_BITS;
            if (sb != 0)
                mpn_rshift(t, Bq + sl, 4 - sl, sb);
            else
                flint_mpn_copyi(t, Bq + sl, 4 - sl);
            s[4] = mpn_add(s, s, 4, t, 4 - sl);
        }
    }
    else
    {
        E = eb;
        sh = 2 * (eb - ea);
        s[4] = 0;
        s[3] = Bq[3]; s[2] = Bq[2]; s[1] = Bq[1]; s[0] = Bq[0];
        if (sh < 4 * FLINT_BITS)
        {
            ulong t[4];
            slong sl = sh / FLINT_BITS, sb = sh % FLINT_BITS;
            if (sb != 0)
                mpn_rshift(t, A + sl, 4 - sl, sb);
            else
                flint_mpn_copyi(t, A + sl, 4 - sl);
            s[4] = mpn_add(s, s, 4, t, 4 - sl);
        }
    }
    derr = (sh > 0);

    /* the normalised three-limb window M 2^X (S >= 2^254: the window
       starts in s4 or s3) */
    if (s[4] != 0)
    {
        m2 = (s[4] << (FLINT_BITS - 1)) | (s[3] >> 1);
        m1 = (s[3] << (FLINT_BITS - 1)) | (s[2] >> 1);
        m0 = (s[2] << (FLINT_BITS - 1)) | (s[1] >> 1);
        X = 2 * E + 1;
    }
    else
    {
        c = flint_clz(s[3]);
        if (c == 0)
        {
            m2 = s[3]; m1 = s[2]; m0 = s[1];
        }
        else
        {
            m2 = (s[3] << c) | (s[2] >> (FLINT_BITS - c));
            m1 = (s[2] << c) | (s[1] >> (FLINT_BITS - c));
            m0 = (s[1] << c) | (s[0] >> (FLINT_BITS - c));
        }
        X = 2 * E - (slong) c;
    }

    if (ARF_SMALL_LOG_OK_2(m2, X))
    {
        re->sgnbit = _arf_small_log_3(&r3, &r2, &r1, &r0, m2, m1, m0, X);
        re->y[3] = r3; re->y[2] = r2; re->y[1] = r1; re->y[0] = r0;
        re->k = -1;
        re->err = ARF_SMALL_LOG_ERR_3 + 6;
        return 0;
    }
    else
    {
        ulong d[4], u[4];
        slong p = 4 * FLINT_BITS - 2 * E, ed, dn;
        int dneg;

        /* D = S - 2^p, p in [3 FLINT_BITS, 4 FLINT_BITS] */
        if (s[4] != 0)
        {
            FLINT_ASSERT(p == 4 * FLINT_BITS);
            d[3] = s[3]; d[2] = s[2]; d[1] = s[1]; d[0] = s[0];
            dneg = 0;
        }
        else if (p == 4 * FLINT_BITS)
        {
            mpn_neg(d, s, 4);
            dneg = 1;
        }
        else
        {
            FLINT_ASSERT(p >= 3 * FLINT_BITS && p < 4 * FLINT_BITS);
            u[3] = UWORD(1) << (p - 3 * FLINT_BITS);
            u[2] = u[1] = u[0] = 0;
            if (mpn_cmp(s, u, 4) >= 0)
            {
                mpn_sub_n(d, s, u, 4);
                dneg = 0;
            }
            else
            {
                mpn_sub_n(d, u, s, 4);
                dneg = 1;
            }
        }

        dn = 4;
        while (dn > 0 && d[dn - 1] == 0)
            dn--;
        if (dn == 0)
            return ARF_ELEM_UNCERTAIN;

        /* the top three limbs of D, normalised (the dropped limb below a
           unit) */
        c = flint_clz(d[dn - 1]);
        {
            ulong t[4];
            if (c != 0)
                mpn_lshift(t, d, dn, c);
            else
                flint_mpn_copyi(t, d, dn);
            m2 = t[dn - 1];
            m1 = (dn >= 2) ? t[dn - 2] : 0;
            m0 = (dn >= 3) ? t[dn - 3] : 0;
            lerr = (dn == 4 && t[0] != 0);
        }
        ed = 2 * E - (4 - dn) * FLINT_BITS - (slong) c;

        if (ed > -24)
            return ARF_ELEM_UNCERTAIN;

        if (derr)
        {
            if (2 * E - ed > 120)
                return ARF_ELEM_UNCERTAIN;
            derr = UWORD(1) << (2 * E - ed);
        }
        derr += lerr;

        _acf_small_log1p_series_3(&r3, &r2, &r1, &r0, &re->err, m2, m1, m0, ed, dneg);
        re->y[3] = r3; re->y[2] = r2; re->y[1] = r1; re->y[0] = r0;
        re->k = ed - 1;
        re->sgnbit = dneg;
        re->err += 2 * derr;
        return 0;
    }
}

/* atan2(b, a) as the frame im (k = 0, or k = -d for a tiny angle) for
   a = (-1)^asgn da 2^(ea - FLINT_BITS), b = (-1)^bsgn db 2^(eb -
   FLINT_BITS): c pi/4 + sigma atan(t) with t = (num / den) 2^-d < 1 the
   smaller of |b| / |a| and |a| / |b| (arf/elem.h); a quotient far below
   1 (d large, where the kernel of arf declines) through one division
   instead, atan t = t - t^3 / 3 + t^5 / 5 - ... with the terms that
   matter */
int
_acf_small_atan2_1(_acf_frame_struct * im, ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn)
{
    ulong num, den, t2, t1, t0, err;
    slong d;
    int c, sigma;

    if (eb < ea || (eb == ea && db <= da))
    {
        num = db; den = da; d = ea - eb;
        c = asgn ? 4 : 0;
        sigma = asgn ? -1 : 1;
    }
    else
    {
        num = da; den = db; d = eb - ea;
        c = 2;
        sigma = asgn ? 1 : -1;
    }

    im->k = 0;
    im->sgnbit = bsgn;

    if (d >= FLINT_BITS || (c == 0 && d > 32))
    {
        ulong q[2], a1, a0;

        /* t = (q / B^2) 2^-d, q = floor(num B^2 / den) with num < den,
           within a unit */
        if (num >= den)
        {
            _mp_real_divq_4_2(q, num >> 1, num << (FLINT_BITS - 1), UWORD(0), UWORD(0), den, UWORD(0));
            d--;
        }
        else
            _mp_real_divq_4_2(q, num, UWORD(0), UWORD(0), UWORD(0), den, UWORD(0));

        if (c == 0)
        {
            /* the frame of t with k = -d: t^3 / 3 is (q / B^2)^3
               2^(2 FLINT_BITS - 2d) / 3 units (a limb, below a unit for
               d >= FLINT_BITS), the tail t^5 / 5 below 2^-6 */
            im->y[2] = 0; im->y[1] = q[1]; im->y[0] = q[0];
            im->k = -d;
            err = 2;
            if (d < FLINT_BITS)
            {
                ulong p1, p0, u;
                _mp_real_sqrhi_2x2(&p1, &p0, q[1], q[0]);
                _mp_real_mulhi_2x2_sloppy(&p1, &p0, p1, p0, q[1], q[0]);
                u = (p1 >> (2 * d - FLINT_BITS)) / 3;
                sub_ddmmss(im->y[1], im->y[0], im->y[1], im->y[0], UWORD(0), u);
                err = 4;
            }
            im->err = err;
            return 0;
        }

        /* t in the units of the result (k = 0): atan t = t within a unit
           for d >= FLINT_BITS */
        _arf_small_rshift_3(&t2, &a1, &a0, UWORD(0), q[1], q[0], d);
        _arf_small_atan_finish_1(&t2, &t1, &t0, &err, a1, a0, 3, c, sigma);
    }
    else
    {
        if (_arf_small_atan2_1(&t2, &t1, &t0, &err, num, den, d, c, sigma) != 0)
            return ARF_ELEM_UNCERTAIN;
    }

    im->y[2] = t2; im->y[1] = t1; im->y[0] = t0;
    im->err = err;
    return 0;
}

int
_acf_small_atan2_2(_acf_frame_struct * im, ulong da1, ulong da0, slong ea, int asgn, ulong db1, ulong db0, slong eb, int bsgn)
{
    ulong num1, num0, den1, den0, t3, t2, t1, t0, err;
    slong d;
    int c, sigma;

    if (eb < ea || (eb == ea && (db1 < da1 || (db1 == da1 && db0 <= da0))))
    {
        num1 = db1; num0 = db0; den1 = da1; den0 = da0; d = ea - eb;
        c = asgn ? 4 : 0;
        sigma = asgn ? -1 : 1;
    }
    else
    {
        num1 = da1; num0 = da0; den1 = db1; den0 = db0; d = eb - ea;
        c = 2;
        sigma = asgn ? 1 : -1;
    }

    im->k = 0;
    im->sgnbit = bsgn;

    if (d >= FLINT_BITS || (c == 0 && d > 40))
    {
        ulong q2, q1, q0, a2, a1, a0;

        if (num1 > den1 || (num1 == den1 && num0 >= den0))
        {
            _mp_real_small_divq_3(&q2, &q1, &q0, num1 >> 1, (num0 >> 1) | (num1 << (FLINT_BITS - 1)),
                num0 << (FLINT_BITS - 1), den1, den0, UWORD(0));
            d--;
        }
        else
            _mp_real_small_divq_3(&q2, &q1, &q0, num1, num0, UWORD(0), den1, den0, UWORD(0));

        if (c == 0)
        {
            /* t^3 / 3 is (q / B^3)^3 2^(3 FLINT_BITS - 2d) / 3 units of the
               frame of t (k = -d), t^5 / 5 is (q / B^3)^5 2^(3 FLINT_BITS
               - 4d) / 5 units: below a unit for d >= 3 FLINT_BITS / 2 resp.
               3 FLINT_BITS / 4; the tail t^7 / 7 below 2^-56 */
            im->y[3] = 0; im->y[2] = q2; im->y[1] = q1; im->y[0] = q0;
            im->k = -d;
            err = 2;
            if (2 * d < 3 * FLINT_BITS)
            {
                ulong p2, p1, p0, w2, w1, w0, u[3];
                _mp_real_small_sqrhi_3x3(&p2, &p1, &p0, q2, q1, q0);
                _mp_real_small_mulhi_3x3(&w2, &w1, &w0, p2, p1, p0, q2, q1, q0);
                _arf_small_rshift_3(u + 2, u + 1, u + 0, w2, w1, w0, 2 * d);
                mpn_divrem_1(u, 0, u, 3, 3);
                sub_dddmmmsss(im->y[2], im->y[1], im->y[0], im->y[2], im->y[1], im->y[0], u[2], u[1], u[0]);
                err += 3;
                if (4 * d < 3 * FLINT_BITS)
                {
                    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, w2, w1, w0, p2, p1, p0);
                    _arf_small_rshift_3(u + 2, u + 1, u + 0, w2, w1, w0, 4 * d);
                    mpn_divrem_1(u, 0, u, 3, 5);
                    add_sssaaaaaa(im->y[2], im->y[1], im->y[0], im->y[2], im->y[1], im->y[0], u[2], u[1], u[0]);
                    err += 3;
                }
            }
            im->err = err;
            return 0;
        }

        _arf_small_rshift_4(&t3, &a2, &a1, &a0, UWORD(0), q2, q1, q0, d);
        _arf_small_atan_finish_2(&t3, &t2, &t1, &t0, &err, a2, a1, a0, 3, c, sigma);
    }
    else
    {
        if (_arf_small_atan2_2(&t3, &t2, &t1, &t0, &err, num1, num0, den1, den0, d, c, sigma) != 0)
            return ARF_ELEM_UNCERTAIN;
    }

    im->y[3] = t3; im->y[2] = t2; im->y[1] = t1; im->y[0] = t0;
    im->err = err;
    return 0;
}

#endif
