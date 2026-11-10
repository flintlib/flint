/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

/*
    sinh, cosh and sinh_cosh (arb/elem.h): the point evaluation
    _arb_sinh_cosh_arf (tiny arguments directly, huge ones by the
    exponential, the register fast paths of arf/elem.h, else the kernel
    _mp_real_sinh_cosh_fp) and the propagation of a radius r through
    |sinh(x + t) - sinh x|, |cosh(x + t) - cosh x| <= |t| cosh(x) e^|t|
    <= r C (1 + 2r) for r <= 1/2 and a bound C on cosh x (from the
    computed cosine, nearly free); wide balls by monotonicity
    (arb_sinh_cosh_wide).
*/

#define MAGLIM(prec) FLINT_MAX(128, 2 * (prec))

static void
arb_sinh_cosh_wide(arb_t s, arb_t c, const arb_t x, slong prec)
{
    mag_t t, u, v, w;

    mag_init(t);
    mag_init(u);
    mag_init(v);
    mag_init(w);

    arb_get_mag_lower(t, x);
    arb_get_mag(u, x);

    if (c != NULL)
    {
        mag_cosh_lower(v, t);
        mag_cosh(w, u);
    }

    if (s != NULL)
    {
        if (mag_is_zero(t))
        {
            arf_get_mag_lower(t, arb_midref(x));
            mag_sub(t, arb_radref(x), t);

            mag_sinh(t, t);
            mag_sinh(u, u);

            if (arf_sgn(arb_midref(x)) > 0)
                arb_set_interval_neg_pos_mag(s, t, u, prec);
            else
                arb_set_interval_neg_pos_mag(s, u, t, prec);
        }
        else
        {
            mag_sinh_lower(t, t);
            mag_sinh(u, u);

            if (arf_sgn(arb_midref(x)) > 0)
            {
                arb_set_interval_mag(s, t, u, prec);
            }
            else
            {
                arb_set_interval_mag(s, t, u, prec);
                arb_neg(s, s);
            }
        }
    }

    if (c != NULL)
        arb_set_interval_mag(c, v, w, prec);

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
    mag_clear(w);
}


/* an upper bound for the frame (y, len) 2^e within err units (err far
   below 2^30) */
static void
_mag_set_frame(mag_t m, nn_srcptr y, slong len, slong e, ulong err)
{
    slong i;
    unsigned int b;
    ulong t;

    for (i = len - 1; i >= 0 && y[i] == 0; i--)
        ;
    if (i < 0)
    {
        mag_set_ui_2exp_si(m, err, e);
        return;
    }
    if (i == 0)
    {
        mag_set_ui_2exp_si(m, y[0], e);
        mag_add_ui_2exp_si(m, m, err + 1, e);
        return;
    }

    /* the top FLINT_BITS bits of (y, i + 1), truncated, then two units of
       its top 30 bits (the truncation and the error) */
    b = FLINT_BIT_COUNT(y[i]);
    t = (b == FLINT_BITS) ? y[i] : ((y[i] << (FLINT_BITS - b)) | (y[i - 1] >> b));
    mag_set_ui_2exp_si(m, (t >> (FLINT_BITS - MAG_BITS)) + 2, e + FLINT_BITS * i + b - MAG_BITS);
}

void
_arb_sinh_cosh_arf(arb_t rs, arb_t rc, arb_t rt, mag_t mc, const arf_t x, slong prec, slong maglim)
{
    nn_srcptr d;
    slong xn, e, nk, w, k;
    ulong err;
    int xsgn;

    if (mc != NULL)
        mag_one(mc);

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            if (rs != NULL) arb_zero(rs);
            if (rc != NULL) arb_one(rc);
            if (rt != NULL) arb_zero(rt);
        }
        else if (arf_is_nan(x))
        {
            if (rs != NULL) arb_indeterminate(rs);
            if (rc != NULL) arb_indeterminate(rc);
            if (rt != NULL) arb_indeterminate(rt);
        }
        else
        {
            int sgn = arf_sgn(x);    /* (the outputs may alias x) */
            if (rs != NULL) arb_set_arf(rs, x);
            if (rc != NULL) arb_pos_inf(rc);
            if (rt != NULL) arb_set_si(rt, sgn);
            if (mc != NULL) mag_inf(mc);
        }
        return;
    }

    xsgn = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e)))
        e = (fmpz_sgn(ARF_EXPREF(x)) < 0) ? -prec - 100 : maglim + 1;

    /* tiny: sinh x = x within |x|^3 / 2, tanh x = x within |x|^3 / 2,
       cosh x = 1 within x^2 (any exponent) */
    if (2 * e <= -prec - 2)
    {
        fmpz_t t3, t2;
        fmpz_init(t3);
        fmpz_init(t2);
        fmpz_mul_ui(t3, ARF_EXPREF(x), 3);
        fmpz_sub_ui(t3, t3, 1);
        fmpz_mul_ui(t2, ARF_EXPREF(x), 2);

        if (rs != NULL || rt != NULL)
        {
            arb_ptr r = (rs != NULL) ? rs : rt;
            int inexact = arf_set_round(arb_midref(r), x, prec, ARB_RND);
            mag_one(arb_radref(r));
            mag_mul_2exp_fmpz(arb_radref(r), arb_radref(r), t3);
            if (inexact)
                arf_mag_add_ulp(arb_radref(r), arb_radref(r), arb_midref(r), prec);
            if (rs != NULL && rt != NULL)
                arb_set(rt, rs);
        }
        if (rc != NULL)
        {
            arf_one(arb_midref(rc));
            mag_one(arb_radref(rc));
            mag_mul_2exp_fmpz(arb_radref(rc), arb_radref(rc), t2);
        }
        if (mc != NULL)
            mag_set_ui_2exp_si(mc, 3, -1);

        fmpz_clear(t3);
        fmpz_clear(t2);
        return;
    }

    if (e > MP_REAL_EXP_FP_MAX_EXP)
    {
        /* sinh |x| and cosh x are exp(|x|) / 2 within a relative
           2 exp(-2|x|), tanh is +-1 within it */
        if (rt != NULL)
        {
            arf_set_si(arb_midref(rt), xsgn ? -1 : 1);
            mag_set_ui_2exp_si(arb_radref(rt), 1, -prec - 10);
        }
        if (rs != NULL || rc != NULL || mc != NULL)
        {
            arb_t t;
            arb_init(t);
            arf_abs(arb_midref(t), x);
            arb_exp_arf(t, arb_midref(t), prec + 2, 0, maglim);
            arb_mul_2exp_si(t, t, -1);
            mag_add_ui_2exp_si(arb_radref(t), arb_radref(t), 1, -prec - 10);
            if (mc != NULL)
                arb_get_mag(mc, t);
            if (rc != NULL)
                arb_set_round(rc, t, prec);
            if (rs != NULL)
            {
                arb_set_round(rs, t, prec);
                if (xsgn)
                    arb_neg(rs, rs);
            }
            arb_clear(t);
        }
        return;
    }

#if FLINT_BITS == 64
    if (xn <= 2 && prec <= 2 * FLINT_BITS && e >= ((xn == 1 && prec <= FLINT_BITS) ? -31 : -40))
    {
        /* (a one-limb mantissa at more than FLINT_BITS bits padded to two) */
        ulong t3, t2, t1, t0, c3, c2, c1, c0;
        ulong d1 = d[xn - 1], d0 = (xn == 2) ? d[0] : 0;
        int one = (xn == 1 && prec <= FLINT_BITS);

        if (rt == NULL)
        {
            if (one)
            {
                _arf_small_sinh_cosh_1(&t2, &t1, &t0, &c2, &c1, &c0, &k, d[0], e);
                if (rs != NULL)
                    _arb_small_set_1(rs, t2, t1, t0, k, ARF_SMALL_HYP_ERR_2, xsgn, prec);
                if (rc != NULL)
                    _arb_small_set_1(rc, c2, c1, c0, k, ARF_SMALL_HYP_ERR_2, 0, prec);
                if (mc != NULL)
                {
                    ulong t[3];
                    t[0] = c0; t[1] = c1; t[2] = c2;
                    _mag_set_frame(mc, t, 3, k - 2 * FLINT_BITS, ARF_SMALL_HYP_ERR_2);
                }
            }
            else
            {
                _arf_small_sinh_cosh_2(&t3, &t2, &t1, &t0, &c3, &c2, &c1, &c0, &k, d1, d0, e);
                if (rs != NULL)
                    _arb_small_set_2(rs, t3, t2, t1, t0, k, ARF_SMALL_HYP_ERR_3, xsgn, prec);
                if (rc != NULL)
                    _arb_small_set_2(rc, c3, c2, c1, c0, k, ARF_SMALL_HYP_ERR_3, 0, prec);
                if (mc != NULL)
                {
                    ulong t[4];
                    t[0] = c0; t[1] = c1; t[2] = c2; t[3] = c3;
                    _mag_set_frame(mc, t, 4, k - 3 * FLINT_BITS, ARF_SMALL_HYP_ERR_3);
                }
            }
            return;
        }
        else
        {
            int code;
            if (one)
            {
                code = _arf_small_tanh_1(&t2, &t1, &t0, &k, &err, d[0], e, ARF_SMALL_HYP_TANH);
                if (code == 0)
                {
                    _arb_small_set_1(rt, t2, t1, t0, k, err, xsgn, prec);
                    return;
                }
            }
            else
            {
                code = _arf_small_tanh_2(&t3, &t2, &t1, &t0, &k, &err, d1, d0, e, ARF_SMALL_HYP_TANH);
                if (code == 0)
                {
                    _arb_small_set_2(rt, t3, t2, t1, t0, k, err, xsgn, prec);
                    return;
                }
            }
            /* (ARF_SMALL_HYP_ONE: +-1 within a unit of the frame) */
            arf_set_si(arb_midref(rt), xsgn ? -1 : 1);
            mag_set_ui_2exp_si(arb_radref(rt), 1, one ? -2 * FLINT_BITS + 4 : -3 * FLINT_BITS + 4);
            return;
        }
    }
#endif

    nk = _arb_elem_limbs(prec);
    {
        nn_ptr ys, yc;
        slong ylen = nk + 2 + ((e < 0) ? (1 - e + FLINT_BITS - 1) / FLINT_BITS : 1);
        TMP_INIT;

        TMP_START;
        ys = TMP_ALLOC(2 * ylen * sizeof(ulong));
        yc = ys + ylen;

        if (rt != NULL)
        {
            _mp_real_tanh_fp(ys, &w, &err, d, xn, e, _arb_elem_bits(prec));
            _arb_set_mpn_err(rt, ys, w + 1, -FLINT_BITS * w, err, err, xsgn, prec);
        }
        else
        {
            k = _mp_real_sinh_cosh_fp((rs != NULL) ? ys : NULL, (rc != NULL || mc != NULL) ? yc : NULL, &w, &err, d, xn, e, _arb_elem_bits(prec));
            if (rs != NULL)
                _arb_set_mpn_err(rs, ys, w + 1, k - FLINT_BITS * w, err, err, xsgn, prec);
            if (rc != NULL)
                _arb_set_mpn_err(rc, yc, w + 1, k - FLINT_BITS * w, err, err, 0, prec);
            if (mc != NULL)
                _mag_set_frame(mc, yc, w + 1, k - FLINT_BITS * w, err);
        }
        TMP_END;
    }
}

/* the propagated error r C (1 + 2r) for r <= 1/2 */
static void
_mag_hyp_propagate(mag_t res, const mag_t C, const mag_t r)
{
    mag_t t;
    mag_init(t);
    mag_mul_2exp_si(t, r, 1);
    mag_add_ui(t, t, 1);
    mag_mul(t, t, C);
    mag_mul(t, t, r);
    mag_add(res, res, t);
    mag_clear(t);
}

/* sinh (s), cosh (c), either NULL */
static void
_arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec), radexp;

    if (mag_is_zero(arb_radref(x)))
    {
        _arb_sinh_cosh_arf(s, c, NULL, NULL, arb_midref(x), prec, maglim);
        return;
    }

    if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
        {
            if (s != NULL) arb_indeterminate(s);
            if (c != NULL) arb_indeterminate(c);
        }
        else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            if (s != NULL) arb_set_arf(s, arb_midref(x));
            if (c != NULL) arb_pos_inf(c);
        }
        else
        {
            if (s != NULL) arb_zero_pm_inf(s);
            if (c != NULL) arb_zero_pm_inf(c);
        }
        return;
    }

    if (mag_cmp_2exp_si(arb_radref(x), -20) > 0)
    {
        if (mag_cmp_2exp_si(arb_radref(x), 10) < 0 && arf_cmpabs_2exp_si(arb_midref(x), 4) < 0)
            arb_sinh_cosh_wide(s, c, x, prec);
        else
        {
            /* a wide ball with a large midpoint: the exponentials */
            arb_t t, u;
            arb_init(t);
            arb_init(u);
            arb_exp_invexp(t, u, x, prec + 4);
            if (s != NULL)
            {
                arb_sub(s, t, u, prec);
                arb_mul_2exp_si(s, s, -1);
            }
            if (c != NULL)
            {
                arb_add(c, t, u, prec);
                arb_mul_2exp_si(c, c, -1);
            }
            arb_clear(t);
            arb_clear(u);
        }
        return;
    }

    /* the precision bounded by the accuracy of the input */
    radexp = MAG_EXP(arb_radref(x));
    if (COEFF_IS_MPZ(radexp))
        radexp = -prec - 100;
    prec = FLINT_MIN(prec, 20 - radexp);
    prec = FLINT_MAX(prec, 2);

    {
        mag_t r, C;
        mag_init_set(r, arb_radref(x));
        mag_init(C);
        _arb_sinh_cosh_arf(s, c, NULL, C, arb_midref(x), prec, maglim);
        if (s != NULL)
            _mag_hyp_propagate(arb_radref(s), C, r);
        if (c != NULL)
            _mag_hyp_propagate(arb_radref(c), C, r);
        mag_clear(r);
        mag_clear(C);
    }
}

void
arb_cosh(arb_t c, const arb_t x, slong prec)
{
    _arb_sinh_cosh(NULL, c, x, prec);
}

void
arb_sinh(arb_t s, const arb_t x, slong prec)
{
    _arb_sinh_cosh(s, NULL, x, prec);
}

void
arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)
{
    if (s == x || c == x)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, x);
        _arb_sinh_cosh(s, c, t, prec);
        arb_clear(t);
    }
    else
        _arb_sinh_cosh(s, c, x, prec);
}
