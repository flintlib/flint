/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "acb_hypgeom.h"

static void
bsplit_zero(acb_t P, acb_t R, acb_t Q, const acb_t z,
    slong a, slong b, slong prec)
{
    if (b - a == 1)
    {
        acb_mul_ui(P, z, a * a, prec);
        acb_set_ui(R, (a + 1) * (a + 1));
        acb_set(Q, R);
    }
    else
    {
        acb_t P2, R2, Q2;
        slong m;

        acb_init(P2);
        acb_init(R2);
        acb_init(Q2);

        m = a + (b - a) / 2;

        bsplit_zero(P, R, Q, z, a, m, prec);
        bsplit_zero(P2, R2, Q2, z, m, b, prec);

        acb_mul(R, R, Q2, prec);
        acb_addmul(R, P, R2, prec);
        acb_mul(P, P, P2, prec);
        acb_mul(Q, Q, Q2, prec);

        acb_clear(P2);
        acb_clear(R2);
        acb_clear(Q2);
    }
}

#if FLINT_BITS == 64
#define DIVLIM 1625              /* just a tuning value */
#define DIVLIM2 1000000000       /* avoid overflow */
#else
#define DIVLIM 40
#define DIVLIM2 30000
#endif

void
acb_hypgeom_dilog_taylor_sum(acb_t res, const acb_t z, slong n, slong prec)
{
    slong k, qk, m, power;
    ulong q;
    acb_t s, t, u;
    acb_ptr zpow;
    int real;

    if (n <= 3)
    {
        if (n <= 1)
            acb_zero(res);
        else if (n == 2)
            acb_set_round(res, z, prec);
        else
        {
            acb_init(t);
            acb_mul(t, z, z, prec);
            acb_mul_2exp_si(t, t, -2);
            acb_add(res, z, t, prec);
            acb_clear(t);
        }
        return;
    }

    /* use binary splitting */
    if (prec > 4000 && acb_bits(z) < prec * 0.02)
    {
        acb_init(s);
        acb_init(t);
        acb_init(u);
        bsplit_zero(s, t, u, z, 1, n, prec);
        acb_add(s, s, t, prec);
        acb_mul(s, s, z, prec);
        acb_div(res, s, u, prec);
        acb_clear(s);
        acb_clear(t);
        acb_clear(u);
        return;
    }

    real = acb_is_real(z);
    k = n - 1;
    m = n_sqrt(n);

    acb_init(s);
    acb_init(t);
    zpow = _acb_vec_init(m + 1);

    _acb_vec_set_powers(zpow, z, m + 1, prec);

    power = (n - 1) % m;

    while (k >= DIVLIM)
    {
        if (k < DIVLIM2)  /* todo: write a div_uiui function? */
        {
            acb_div_ui(t, zpow + power, k * k, prec);
        }
        else
        {
            acb_div_ui(t, zpow + power, k, prec);
            acb_div_ui(t, t, k, prec);
        }

        acb_add(s, s, t, prec);

        if (power == 0)
        {
            acb_mul(s, s, zpow + m, prec);
            power = m - 1;
        }
        else
        {
            power--;
        }

        k--;
    }

    qk = k;   /* k at which to change denominator */
    q = 1;

    while (k >= 2)
    {
        /* find next qk such that the consecutive denominators can be
           collected in a single word */
        if (k == qk)
        {
            if (q != 1)
                acb_div_ui(s, s, q, prec);

            q = qk * qk;
            qk--;


            while (qk > 1)
            {
                ulong hi, lo;
                umul_ppmm(hi, lo, q, qk * qk);
                if (hi != 0)
                    break;
                q *= qk * qk;
                qk--;
            }

            acb_mul_ui(s, s, q, prec);
        }

        if (power == 0)
        {
            acb_add_ui(s, s, q / (k * k), prec);
            acb_mul(s, s, zpow + m, prec);
            power = m - 1;
        }
        else
        {
            if (real)  /* minor optimization */
                arb_addmul_ui(acb_realref(s), acb_realref(zpow + power), q / (k * k), prec);
            else
                acb_addmul_ui(s, zpow + power, q / (k * k), prec);
            power--;
        }

        k--;
    }

    if (q != 1)
        acb_div_ui(s, s, q, prec);
    acb_add(s, s, z, prec);
    acb_swap(res, s);

    _acb_vec_clear(zpow, m + 1);
    acb_clear(s);
    acb_clear(t);
}

void
acb_hypgeom_dilog_zero_taylor(acb_t res, const acb_t z, slong prec)
{
    mag_t m;
    slong n;
    double x;
    int real;

    mag_init(m);
    acb_get_mag(m, z);
    real = acb_is_real(z);
    x = -mag_get_d_log2_approx(m);

    n = 2;
    if (x > 0.01)
    {
        n = prec / x + 1;
        n += (x > 2.0);  /* relative error for very small |z| */
    }

    n = FLINT_MAX(n, 2);

    mag_geom_series(m, m, n);
    mag_div_ui(m, m, n);
    mag_div_ui(m, m, n);

    if (mag_is_finite(m))
    {
        acb_hypgeom_dilog_taylor_sum(res, z, n, prec);
        if (real)
            arb_add_error_mag(acb_realref(res), m);
        else
            acb_add_error_mag(res, m);
    }
    else
    {
        acb_indeterminate(res);
    }

    mag_clear(m);
}

