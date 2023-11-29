/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "thread_support.h"
#include "arb.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

/* Computes sin(x) or cos(x) using Taylor series truncated at x^N exclusive.
   Computes error bound automatically. Does not allow aliasing of s and x.  */
void
arb_sin_cos_taylor_sum_rs(arb_t s, const arb_t x, slong N, int cosine, slong prec)
{
    mag_t err;
    mag_init(err);

    arb_get_mag(err, x);
    mag_exp_tail(err, err, N);

    if (N == 0 || (!cosine && N == 1))
    {
        arb_zero(s);
    }
    else if (cosine && N <= 2)
    {
        arb_one(s);
    }
    else if (!cosine && N <= 3) /* x */
    {
        arb_set_round(s, x, prec);
    }
    else if (cosine && N <= 4)  /* 1 - x^2/2 */
    {
        arb_mul(s, x, x, prec / 2 + 4);
        arb_mul_2exp_si(s, s, -1);
        arb_sub_ui(s, s, 1, prec);
        arb_neg(s, s);
    }
    else if (!cosine && N <= 5) /* x - x^3/6 */
    {
        arb_mul(s, x, x, prec / 2 + 4);
        arb_div_ui(s, s, 6, prec / 2 + 4);
        arb_mul(s, s, x, prec / 2 + 4);
        arb_sub(s, x, s, prec);
    }
    else
    {
        arb_ptr tpow;
        slong j, k, m, M, tp, xmag;
        mp_limb_t c, d, chi, clo;

        xmag = arf_abs_bound_lt_2exp_si(arb_midref(x));

        /* Convert to order as a series in x^2. */
        if (cosine)
            M = (N + 1) / 2;
        else
            M = N / 2;

        m = n_sqrt(M);

        /* not intended (and not 32-bit safe...) */
        if (M > 30000)
        {
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
        }

        tpow = _arb_vec_init(m + 2);

        arb_mul(s, x, x, prec);
        _arb_vec_set_powers(tpow, s, m + 1, prec);
        arb_zero(s);

        c = 1;
        j = (M - 1) % m;

        for (k = M - 1; k >= 0; k--)
        {
            tp = prec - 2 * k * (-xmag) + 10;
            tp = FLINT_MAX(tp, 2);
            tp = FLINT_MIN(tp, prec);

            if (cosine)
                d = (2 * k - 1) * (2 * k);
            else
                d = (2 * k) * (2 * k + 1);

            if (k != 0)
            {
                umul_ppmm(chi, clo, c, d);

                if (chi != 0)
                {
                    arb_div_ui(s, s, c, tp);
                    c = 1;
                }
            }

            if (k % 2 == 0)
                arb_addmul_ui(s, tpow + j, c, tp);
            else
                arb_submul_ui(s, tpow + j, c, tp);

            if (k != 0)
            {
                c *= d;

                if (j == 0)
                {
                    if (tp > 300000)
                    {
                        arb_set_round(tpow + m + 1, tpow + m, tp);
                        arb_mul(s, s, tpow + m + 1, tp);
                    }
                    else
                    {
                        arb_mul(s, s, tpow + m, tp);
                    }

                    j = m - 1;
                }
                else
                {
                    j--;
                }
            }
        }

        arb_div_ui(s, s, c, prec);
        if (!cosine)
            arb_mul(s, s, x, prec);

        _arb_vec_clear(tpow, m + 2);
    }

    arb_add_error_mag(s, err);
    mag_clear(err);
}

void
arb_sin_cos_arf_rs_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec)
{
    slong q, xmag, wp, k, N;
    arb_t s, t;
    int negate;

    if (arf_is_zero(x))
    {
        if (res_sin != NULL)
            arb_zero(res_sin);
        if (res_cos != NULL)
            arb_one(res_cos);
        return;
    }

    xmag = arf_abs_bound_lt_2exp_si(x);

    /* x + O(x^3), 1 + O(x^2) */
    if (xmag < -(prec / 2) - 4)
    {
        arb_init(t);
        arf_set(arb_midref(t), x);
        if (res_sin != NULL)
            arb_sin_cos_taylor_sum_rs(res_sin, t, 3, 0, prec);
        if (res_cos != NULL)
            arb_sin_cos_taylor_sum_rs(res_cos, t, 2, 1, prec);
        arb_clear(t);
        return;
    }

    xmag = FLINT_MAX(xmag, -prec);

    /* could include sanity test, but we assume that the function
       gets called with a valid |x| < pi/2
    if (xmag > 0)
    {
        if (res_sin != NULL)
            arb_indeterminate(res_sin);
        if (res_cos != NULL)
            arb_indeterminate(res_cos);
        return;
    }
     */

    arb_init(s);
    arb_init(t);

    negate = arf_sgn(x) < 0;

    /* generic tuning value */
    q = 4.5 * pow(prec, 0.2);
    q = FLINT_MAX(q, 6);
    /* adjust to magnitude */
    q = FLINT_MAX(0, xmag + q);
    /* don't do a redundant square root */
    if (q <= 2)
        q = 0;

    wp = prec + 10 + 2 * FLINT_BIT_COUNT(prec);

    /* t = x/2^q */
    arf_mul_2exp_si(arb_midref(t), x, -q);

    if (q == 0 && res_sin != NULL)
    {
        /* compute cos from sin since the square root has less cancellation */
        wp += (-xmag);
        N = _arb_exp_taylor_bound(xmag, wp);
        arb_sin_cos_taylor_sum_rs(s, t, N, 0, wp);

        if (res_sin != NULL)
            arb_set_round(res_sin, s, prec);

        if (res_cos != NULL)
        {
            arb_mul(t, s, s, wp);
            arb_sub_ui(t, t, 1, wp);
            arb_neg(t, t);
            arb_sqrt(res_cos, t, prec);
        }
    }
    else
    {
        /* compute sin from cos */
        wp = prec + 10 + 2 * FLINT_BIT_COUNT(prec);
        wp += 2 * (q - xmag);  /* todo: too much when only computing cos? */

        N = _arb_exp_taylor_bound(xmag - q, wp);

        arb_sin_cos_taylor_sum_rs(s, t, N, 1, wp);

        for (k = 0; k < q; k++)
        {
            arb_mul(s, s, s, wp);
            arb_mul_2exp_si(s, s, 1);
            arb_sub_ui(s, s, 1, wp);
        }

        if (res_cos != NULL)
            arb_set_round(res_cos, s, prec);

        if (res_sin != NULL)
        {
            arb_mul(s, s, s, wp);
            arb_sub_ui(s, s, 1, wp);
            arb_neg(s, s);
            arb_sqrtpos(res_sin, s, prec);
            if (negate)
                arb_neg(res_sin, res_sin);
        }
    }

    arb_clear(s);
    arb_clear(t);
}

void
arb_sin_cos_arf_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec)
{
    arb_t pi4, t, u, v;
    fmpz_t q;
    slong wp, mag;
    int octant, swapsincos, sinnegative, cosnegative, negative;

    mag = arf_abs_bound_lt_2exp_si(x);

    if (mag > FLINT_MAX(65536, 4 * prec))
    {
        if (res_sin != NULL)
            arb_zero_pm_one(res_sin);
        if (res_cos != NULL)
            arb_zero_pm_one(res_cos);
    }
    else if (mag <= 0)  /* todo: compare with pi/4-eps instead? */
    {
        int want_rs;

        if (prec < 20000 || mag < -prec / 16)
        {
            want_rs = 1;
        }
        else if (arf_bits(x) < prec / 128)
        {
            want_rs = 0;
        }
        else if (flint_get_num_available_threads() == 1)
        {
            want_rs = (prec < 200000) || (prec < 1000000000 && mag < -prec / 5000);
        }
        else
        {
            want_rs = (prec < 20000) || (prec < 1000000000 && mag < -prec / 400);
        }

        if (want_rs)
            arb_sin_cos_arf_rs_generic(res_sin, res_cos, x, prec);
        else
            arb_sin_cos_arf_bb(res_sin, res_cos, x, prec);
    }
    else
    {
        arb_init(pi4);
        arb_init(t);
        arb_init(u);
        arb_init(v);
        fmpz_init(q);

        wp = prec + mag + 10;
        negative = arf_sgn(x) < 0;

        arb_const_pi(pi4, wp);
        arb_mul_2exp_si(pi4, pi4, -2);
        arb_set_arf(t, x);
        arb_abs(t, t);

        arb_set_round(v, t, mag + 10);
        arb_set_round(u, pi4, mag + 10);
        arb_div(u, v, u, mag + 10);

        arf_get_fmpz(q, arb_midref(u), ARF_RND_DOWN);
        arb_submul_fmpz(t, pi4, q, wp);
        octant = fmpz_fdiv_ui(q, 8);
        if (octant & 1)
            arb_sub(t, pi4, t, wp);

        arb_clear(pi4);
        arb_clear(u);
        arb_clear(v);

        sinnegative = (octant >= 4) ^ negative;
        cosnegative = (octant >= 2 && octant <= 5);
        swapsincos = (octant == 1 || octant == 2 || octant == 5 || octant == 6);

        /* guard against infinite recursion */
        if (arf_cmpabs_2exp_si(arb_midref(t), 0) > 0)
        {
            flint_throw(FLINT_ERROR, "mod pi/4 reduction unexpectedly failed!\n");
        }

        /* todo: allow NULL in arb_sin_cos and simplify here */
        if (swapsincos)
        {
            if (res_sin != NULL && res_cos != NULL)
                arb_sin_cos(res_cos, res_sin, t, prec);
            else if (res_sin != NULL)
                arb_cos(res_sin, t, prec);
            else
                arb_sin(res_cos, t, prec);
        }
        else
        {
            if (res_sin != NULL && res_cos != NULL)
                arb_sin_cos(res_sin, res_cos, t, prec);
            else if (res_sin != NULL)
                arb_sin(res_sin, t, prec);
            else
                arb_cos(res_cos, t, prec);
        }

        if (sinnegative && res_sin != NULL)
            arb_neg(res_sin, res_sin);
        if (cosnegative && res_cos != NULL)
            arb_neg(res_cos, res_cos);

        arb_clear(t);
        fmpz_clear(q);
    }
}

