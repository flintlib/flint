/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

static void
evaluate_rect(acb_t res, const short * term_prec, slong len, const acb_t x, slong prec)
{
    slong i, j, m, r, n1, n2;
    acb_ptr xs;
    acb_t s, t;
    arb_struct c[17];

    m = n_sqrt(len) + 1;
    m = FLINT_MIN(m, 16);
    r = (len + m - 1) / m;

    xs = _acb_vec_init(m + 1);
    acb_init(s);
    acb_init(t);
    _acb_vec_set_powers(xs, x, m + 1, prec);

    acb_zero(res);

    for (i = r - 1; i >= 0; i--)
    {
        n1 = m * i;
        n2 = FLINT_MIN(len, n1 + m);

        for (j = n1; j < n2; j++)
        {
            if (j == 0)
            {
                arb_init(c);
                arb_one(c);
            }
            else
            {
                if (!_arb_hypgeom_gamma_coeff_shallow(arb_midref(c + j - n1), arb_radref(c + j - n1), j, term_prec[j]))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
            }
        }

        arb_dot(acb_realref(s), NULL, 0, acb_realref(xs), 2, c, 1, n2 - n1, prec);
        arb_dot(acb_imagref(s), NULL, 0, acb_imagref(xs), 2, c, 1, n2 - n1, prec);

#if 0
        acb_set_round(t, xs + m, term_prec[n1]);
        acb_mul(res, res, t, term_prec[n1]);
        acb_add(res, res, s, term_prec[n1]);
#else
        acb_mul(res, res, xs + m, term_prec[n1]);
        acb_add(res, res, s, term_prec[n1]);
#endif
    }

    _acb_vec_clear(xs, m + 1);
    acb_clear(s);
    acb_clear(t);
}

/* Bound requires: |u| <= 20, N <= 10000, N != (1443, 2005, 9891). */
static void
error_bound(mag_t err, const acb_t u, slong N)
{
    mag_t t;
    mag_init(t);

    acb_get_mag(t, u);

    if (N >= 1443 || mag_cmp_2exp_si(t, 4) > 0)
    {
        mag_inf(err);
    }
    else
    {
        mag_pow_ui(err, t, N);
        mag_mul_2exp_si(err, err, arb_hypgeom_gamma_coeffs[N].exp);

        if (mag_cmp_2exp_si(t, -1) > 0)
            mag_mul(err, err, t);
        else
            mag_mul_2exp_si(err, err, -1);

        mag_mul_2exp_si(err, err, 3);

        if (mag_cmp_2exp_si(err, -8) > 0)
            mag_inf(err);
    }

    mag_clear(t);
}

static double
want_taylor(double x, double y, slong prec)
{
    if (y < 0.0) y = -y;
    if (x < 0.0) x = -x;

    if ((prec < 128 && y > 4.0) || (prec < 256 && y > 5.0) ||
        (prec < 512 && y > 8.0) || (prec < 1024 && y > 9.0) || y > 10.0)
    {
        return 0;
    }

    if (x * (1.0 + 0.75 * y) > 8 + 0.15 * prec)
    {
        return 0;
    }

    return 1;
}

int
acb_hypgeom_gamma_taylor(acb_t res, const acb_t z, int reciprocal, slong prec)
{
    acb_t s, u;
    int success;
    double dua, dub, du2, log2u;
    slong i, r, n, wp, tail_bound, goal;
    short term_prec[ARB_HYPGEOM_GAMMA_TAB_NUM];
    mag_t err;

    if (!acb_is_finite(z) ||
        arf_cmp_2exp_si(arb_midref(acb_imagref(z)), 4) >= 0 ||
        arf_cmp_2exp_si(arb_midref(acb_realref(z)), 10) >= 0)
    {
        return 0;
    }

    dua = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_UP);
    dub = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_UP);
    dub = fabs(dub);

    if (!want_taylor(dua, dub, prec))
        return 0;

    if (dua >= 0.0)
        r = (slong) (dua + 0.5);
    else
        r = -(slong) (-dua + 0.5);

    acb_init(s);
    acb_init(u);
    mag_init(err);

    success = 0;

    /* Argument reduction: u = z - r */
    acb_sub_si(u, z, r, 2 * prec + 10);
    dua -= r;

    goal = acb_rel_accuracy_bits(u);

    /* not designed for wide intervals (yet) */
    if (goal < 8)
    {
        success = 0;
        goto cleanup;
    }

    goal = FLINT_MIN(goal, prec - MAG_BITS) + MAG_BITS;
    goal = FLINT_MAX(goal, 5);
    goal = goal + 5;
    wp = goal + 4 + FLINT_BIT_COUNT(FLINT_ABS(r));

    if (wp > ARB_HYPGEOM_GAMMA_TAB_PREC)
    {
        success = 0;
        goto cleanup;
    }

    if (!want_taylor(r, dub, goal))
    {
        success = 0;
        goto cleanup;
    }

    du2 = dua * dua + dub * dub;

    if (du2 > 1e-8)
    {
        log2u = 0.5 * mag_d_log_upper_bound(du2) * 1.4426950408889634074 * (1 + 1e-14);
    }
    else
    {
        slong aexp, bexp;

        aexp = arf_cmpabs_2exp_si(arb_midref(acb_realref(u)), -wp) >= 0 ? ARF_EXP(arb_midref(acb_realref(u))) : -wp;
        bexp = arf_cmpabs_2exp_si(arb_midref(acb_imagref(u)), -wp) >= 0 ? ARF_EXP(arb_midref(acb_imagref(u))) : -wp;
        log2u = FLINT_MAX(aexp, bexp) + 1;
    }

    term_prec[0] = wp;
    n = 0;

    for (i = 1; i < ARB_HYPGEOM_GAMMA_TAB_NUM; i++)
    {
        tail_bound = arb_hypgeom_gamma_coeffs[i].exp + i * log2u + 5;

        if (tail_bound <= -goal)
        {
            n = i;
            break;
        }

        term_prec[i] = FLINT_MIN(FLINT_MAX(wp + tail_bound, 2), wp);

        if (term_prec[i] > arb_hypgeom_gamma_coeffs[i].nlimbs * FLINT_BITS)
        {
            success = 0;
            goto cleanup;
        }
    }

    if (n != 0)
        error_bound(err, u, n);

    if (n == 0 || mag_is_inf(err))
    {
        success = 0;
        goto cleanup;
    }

    evaluate_rect(s, term_prec, n, u, wp);
    acb_add_error_mag(s, err);

    if (r == 0 || r == 1)
    {
        if (r == 0)
            acb_mul(s, s, u, wp);

        if (reciprocal)
        {
            acb_set_round(res, s, prec);
        }
        else
        {
            acb_one(u);
            acb_div(res, u, s, prec);
        }
    }
    else if (r >= 2)
    {
        acb_add_ui(u, u, 1, wp);
        acb_hypgeom_rising_ui_rec(u, u, r - 1, wp);

        if (reciprocal)
            acb_div(res, s, u, prec);
        else
            acb_div(res, u, s, prec);
    }
    else
    {
        /* gamma(x) = (-1)^r / (rgamma(1+x-r)*rf(1+r-x,-r)*(x-r)) */
        /* 1/gamma(x) = (-1)^r * rgamma(1+x-r) * rf(1+r-x,-r) * (x-r) */

        acb_neg(res, z);
        acb_add_si(res, res, 1 + r, wp);
        acb_hypgeom_rising_ui_rec(res, res, -r, wp);
        acb_mul(u, res, u, wp);

        if (reciprocal)
        {
            acb_mul(res, s, u, prec);
        }
        else
        {
            acb_mul(u, s, u, wp);
            acb_inv(res, u, prec);
        }

        if (r % 2)
            acb_neg(res, res);
    }

    success = 1;

cleanup:
    acb_clear(s);
    acb_clear(u);
    mag_clear(err);

    return success;
}

