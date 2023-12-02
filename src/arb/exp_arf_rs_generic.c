/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
arb_exp_taylor_sum_rs_generic(arb_t res, const arb_t x, slong N, slong prec)
{
    arb_t s;
    mag_t err;

    arb_init(s);
    mag_init(err);

    arb_get_mag(err, x);
    mag_exp_tail(err, err, N);

    if (N <= 2)
    {
        if (N == 0)
            arb_zero(res);
        else if (N == 1)
            arb_one(res);
        else if (N == 2)
            arb_add_ui(res, x, 1, prec);
        arb_add_error_mag(res, err);
    }
    else
    {
        arb_ptr tpow;
        slong j, k, m, M, tp, xmag;
        mp_limb_t c, d, chi, clo;

        xmag = arf_abs_bound_lt_2exp_si(arb_midref(x));

        /* Convert to order as a series in x^2. */
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

            arb_addmul_ui(s, tpow + j, c, tp);

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
        arb_mul(s, s, x, prec);

        arb_add_error_mag(s, err);

        /* exp = sinh + sqrt(1 + sinh^2) */
        arb_mul(res, s, s, prec);
        arb_add_ui(res, res, 1, prec);
        arb_sqrt(res, res, prec);
        arb_add(res, res, s, prec);

        _arb_vec_clear(tpow, m + 2);
    }

    mag_clear(err);
    arb_clear(s);
}

void
arb_exp_arf_rs_generic(arb_t res, const arf_t x, slong prec, int minus_one)
{
    slong q, xmag, wp, k, N;
    arb_t t;

    if (arf_is_zero(x))
    {
        if (minus_one)
            arb_zero(res);
        else
            arb_one(res);
        return;
    }

    xmag = arf_abs_bound_lt_2exp_si(x);

    /* 1 + x + O(x^2) */
    /* We don't really need to worry too much about degenerate input
       because the main exp function already takes care of it. */
    if (xmag < -prec - 4)
    {
        mag_t t;
        mag_init(t);
        arf_get_mag(t, x);
        mag_exp_tail(t, t, 2);
        arb_set_arf(res, x);
        arb_add_ui(res, res, minus_one ? 0 : 1, prec);
        arb_add_error_mag(res, t);
        mag_clear(t);
        return;
    }

    arb_init(t);

    /* generic tuning value */
    q = 4.5 * pow(prec, 0.2);
    q = FLINT_MAX(q, 6);
    /* adjust to magnitude */
    q = FLINT_MAX(0, xmag + q);

    wp = prec + 10 + 2 * q + 2 * FLINT_BIT_COUNT(prec);
    if (minus_one && xmag < 0)
        wp += (-xmag);

    /* t = x/2^q */
    arf_mul_2exp_si(arb_midref(t), x, -q);

    N = _arb_exp_taylor_bound(xmag - q, wp);
    arb_exp_taylor_sum_rs_generic(t, t, N, wp);

    /* exp(x) = exp(x/2^q)^(2^q) */
    for (k = 0; k < q; k++)
        arb_mul(t, t, t, wp);

    if (minus_one)
        arb_sub_ui(t, t, 1, wp);

    arb_set_round(res, t, prec);
    arb_clear(t);
}

