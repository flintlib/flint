/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
# define fabs __builtin_fabs
# define log __builtin_log
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

/* todo: improve for small k */
static double log2_bin_uiui_fast(ulong n, ulong k)
{
    static const float htab[] = {0.2007, 0.3374, 0.4490, 0.5437, 0.6254, 0.6963,
        0.7580, 0.8114, 0.8572, 0.8961, 0.9285, 0.9545, 0.9746,
        0.9888, 0.9973, 1.0};

    if (k == 0 || k >= n)
        return 0;
    if (k > n / 2)
        k = n - k;
    k = (32.0 * k) / n;
    return n * htab[FLINT_MIN(k, 15)];
}

/* x is unused but part of the API */
void
arb_hypgeom_legendre_p_ui_deriv_bound(mag_t dp, mag_t dp2, ulong n, const arb_t x, const arb_t x2sub1)
{
    mag_t c, t, u;

    mag_init(c);
    mag_init(t);
    mag_init(u);

    arb_get_mag_lower(t, x2sub1);
    mag_rsqrt(t, t);                /* t >= 1/(1-x^2)^(1/2) */
    mag_mul_ui(u, t, n);            /* u >= n/(1-x^2)^(1/2) */
    mag_set_ui_2exp_si(c, 409, -8); /* c >= 2^(3/2)/sqrt(pi) */
    mag_sqrt(dp, u);
    mag_mul(dp, dp, t);
    mag_mul(dp, dp, c);             /* dp >= c*sqrt(n)/(1-x^2)^(3/4) */

    mag_mul(dp2, dp, u);
    mag_mul_2exp_si(dp2, dp2, 1);   /* dp2 >= 2*c*n^(3/2)/(1-x^2)^(5/4) */

    mag_set_ui(t, n);
    mag_add_ui(t, t, 1);
    mag_mul(t, t, t);               /* t >= (n+1)^2 */
    mag_mul_2exp_si(u, t, -1);      /* u >= (n+1)^2/2 */
    mag_min(dp, dp, u);             /* |P'(x)| <= dp */

    mag_mul(t, t, t);
    mag_mul_2exp_si(u, t, -3);
    mag_min(dp2, dp2, u);           /* |P''(x)| <= dp2 */

    mag_clear(c);
    mag_clear(t);
    mag_clear(u);
}

void
arb_hypgeom_legendre_p_ui(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)
{
    arb_t xsub1, x2sub1;
    double xx, xxsub1, cancellation_zero, cancellation_one;
    double cost_zero, cost_one, cost_asymp;
    double log2x, log2u, tolerance, asymp_error;
    double yy, log2nsy, log2k, size;
    slong wp;
    slong d, k, K_zero, K_one, K_asymp;
    int basecase_ok;

    if (!arb_is_finite(x) || n > UWORD_MAX / 4)
    {
        if (res != NULL)
            arb_indeterminate(res);
        if (res_prime != NULL)
            arb_indeterminate(res_prime);
        return;
    }

    if (arf_sgn(arb_midref(x)) < 0)
    {
        arb_t t;
        arb_init(t);
        arb_neg(t, x);
        arb_hypgeom_legendre_p_ui(res, res_prime, n, t, prec);
        if (n % 2 == 1 && res != NULL)
            arb_neg(res, res);
        if (n % 2 == 0 && res_prime != NULL)
            arb_neg(res_prime, res_prime);
        arb_clear(t);
        return;
    }

    if (arb_is_one(x) && n < UWORD_MAX)
    {
        if (res != NULL)
            arb_set(res, x);
        if (res_prime != NULL)
        {
            arb_set_ui(res_prime, n);
            arb_mul_ui(res_prime, res_prime, n + 1, prec);
            arb_mul_2exp_si(res_prime, res_prime, -1);
        }
        return;
    }

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    if (n == 1)
    {
        if (res != NULL) arb_set_round(res, x, prec);
        if (res_prime != NULL) arb_one(res_prime);
        return;
    }

    xx = arf_get_d(arb_midref(x), ARF_RND_UP);

    /* Use basecase recurrence? */
    /* The following tests are not very elegant, and not completely accurate
       either, but they are fast in the common case. */
    if (res_prime != NULL)
    {
        basecase_ok = ((xx < 0.999999 && n < 10 && prec < 2000) ||
                       (xx < 0.999999 && n < 50 && prec < 1000) ||
                       (xx < 0.9999 && n < 100 && prec < 1000) ||
                       (xx < 0.999  && n < 350 && prec < 1000) ||
                       (xx < 0.9 && n < 400 && prec < 1000))
                   && ((xx > 0.00001 && n < 10 && prec < 2000) ||
                       (xx > 0.00001 && n < 60 && prec < 1000) ||
                       (xx > 0.01 && n < 200 && prec < 1000) ||
                       (xx > 0.1 && n < 400 && prec < 1000));

        /* the recurrence also performs better when n ~= prec */
        if (!basecase_ok)
            basecase_ok = (xx > 0.1 && xx < 0.99 && n < 800
                && prec > 0.4 * n && prec < 1.5 * n);
    }
    else if (prec < 500)
    {
        basecase_ok = ((xx < 0.999999 && n < 20) ||
                       (xx < 0.999 && n < 60) ||
                       (xx < 0.9 && n < 100))
                   && ((xx > 0.00001 && n < 20) ||
                       (xx > 0.01 && n < 60) ||
                       (xx > 0.1 && n < 100));

        if (!basecase_ok)
            basecase_ok = (xx > 0.1 && xx < 0.99 && n < 300
                && prec > 0.4 * n && prec < 1.5 * n);
    }
    else
    {
        basecase_ok = 0;
    }

    if (basecase_ok)
    {
        mag_t t;
        mag_init(t);
        arb_get_mag(t, x);
        if (mag_cmp_2exp_si(t, 0) >= 0)
            basecase_ok = 0;
        mag_clear(t);
    }

    if (basecase_ok)
    {
        arb_hypgeom_legendre_p_ui_rec(res, res_prime, n, x, prec);
        return;
    }

    arb_init(xsub1);
    arb_init(x2sub1);

    arb_sub_ui(xsub1, x, 1, prec + 10);

    arb_mul(x2sub1, x, x, 2 * prec);
    arb_sub_ui(x2sub1, x2sub1, 1, prec + 10);
    arb_neg(x2sub1, x2sub1);

    /* use series at 1 unless |x| < 1-eps */
    if (!arb_is_negative(xsub1) ||
        arf_cmp_d(arb_midref(xsub1), ldexp(1.0, -2 * FLINT_BIT_COUNT(n))) > 0)
    {
        if (arf_cmp_d(arb_midref(xsub1), 2.0) >= 0)
        {
            if (n < 10000.0 * prec && n < UWORD_MAX / 4)
                K_one = n + 1;
            else
                K_one = 1;
        }
        else  /* check for early convergence */
        {
            xxsub1 = arf_get_d(arb_midref(xsub1), ARF_RND_UP);
            log2u = log(fabs(xxsub1) * 0.5) * 1.44269504088896;
            if (log2u < -30)
                log2u = arf_abs_bound_lt_2exp_si(arb_midref(xsub1)) - 1.0;

            K_one = n + 1;
            K_one = FLINT_MIN(K_one, 100000.0 * prec);
            K_one = FLINT_MIN(K_one, UWORD_MAX / 4);

            size = 0.0;

            if (n * (2.0 + log2u) < -prec)
            {
                for (k = 1; k < K_one; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, k)
                        + log2_bin_uiui_fast(n + k, k) + k * log2u;

                    if (size < -prec)
                    {
                        K_one = k;
                        break;
                    }
                }
            }
        }

        arb_hypgeom_legendre_p_ui_one(res, res_prime, n, x, K_one, prec);
    }
    else   /* guaranteed to have |x| < 1 */
    {
        cost_zero = 1e100;
        cost_one = 1e100;
        cost_asymp = 1e100;

        xx = FLINT_MAX(xx, 1e-50);
        xxsub1 = arf_get_d(arb_midref(xsub1), ARF_RND_UP);

        /* Estimate cancellation for series expansion at 0. */
        /* |P_n(xi)| ~= (x+sqrt(1+x^2))^n. */
        cancellation_zero = n * log(xx + sqrt(1.0 + xx * xx)) * 1.44269504088896;
        cancellation_zero = FLINT_MIN(cancellation_zero, 1.272 * n);
        cancellation_zero = FLINT_MAX(cancellation_zero, 0.0);

        /* Estimate cancellation for series expansion at 1. */
        /* For x >= 1, P_n(x) ~= I_0(n*sqrt(2(x-1))) ~= exp(n*sqrt(2(x-1))) */
        if (xxsub1 >= 0.0)
        {
            cancellation_one = 0.0;
        }
        else
        {
            cancellation_one = n * sqrt(2.0*fabs(xxsub1)) * 1.44269504088896;
            cancellation_one = FLINT_MIN(cancellation_one, 2.0 * n);
            cancellation_one = FLINT_MAX(cancellation_one, 0.0);
        }

        d = n / 2;
        K_zero = d + 1;
        K_one = n + 1;
        K_asymp = 1;
        asymp_error = 0.0;

        wp = 1.01 * prec + FLINT_BIT_COUNT(n);
        tolerance = -wp;

        /* adjust for relative tolerance near 0 */
        if (n % 2)
        {
            tolerance += arf_abs_bound_lt_2exp_si(arb_midref(x));
        }

        if (n > 10)
        {
            /* look for early truncation of series at 1 */
            log2u = log(fabs(xxsub1) * 0.5) * 1.44269504088896;
            if (log2u < -30)
                log2u = arf_abs_bound_lt_2exp_si(arb_midref(xsub1)) - 1.0;

            log2x = log(fabs(xx)) * 1.44269504088896;
            if (log2x < -30)
                log2x = arf_abs_bound_lt_2exp_si(arb_midref(x));

            if (n * (2.0 + log2u) < tolerance)
            {
                for (k = 1; k < K_one; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, k)
                        + log2_bin_uiui_fast(n + k, k) + k * log2u;

                    if (size < tolerance)
                    {
                        K_one = k;
                        break;
                    }
                }
            }

            /* look for early truncation of series at 0 */
            if (n * (1.0 + log2x) < tolerance)
            {
                for (k = 1; k < K_zero; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, d - k)
                        + log2_bin_uiui_fast(n+1+2*k, n) - n + 2.0 * k * log2x;

                    if (size < tolerance)
                    {
                        K_zero = k;
                        break;
                    }
                }
            }

            /* look for possible convergence of asymptotic series */
            /* requires information about y = sqrt(1-x^2) */
            yy = arf_get_d(arb_midref(x2sub1), ARF_RND_DOWN);
            yy = sqrt(FLINT_MAX(0.0, yy));
            log2nsy = log(2.0 * n * yy) * 1.44269504088896;

            cost_zero = (prec + cancellation_zero) * K_zero;
            cost_one = (prec + cancellation_one) * K_one;

            for (k = 1; k < n &&
                    prec * k < FLINT_MIN(cost_zero, cost_one);
                        k = FLINT_MAX(k + 1, k * 1.05))
            {
                /* todo: better account for prefactor in the asymptotic series? */
                log2k = log(k) * 1.44269504088896;
                size = 3.0 + k * (log2k - 1.43);  /* estimate K! */
                size -= k * log2nsy;              /* 1/(2n sin(theta))^K */

                if (size < asymp_error)
                {
                    asymp_error = size;
                    K_asymp = k;
                }

                if (size < tolerance)
                {
                    break;
                }
            }
        }

        cost_zero = (prec + cancellation_zero) * K_zero;
        cost_one = (prec + cancellation_one) * K_one;
        cost_asymp = (prec + 0.0) * K_asymp * 2.0;

#if 0
        printf("zero:  K = %ld, cost = %g, cancel %g\n",                K_zero, cost_zero, cancellation_zero);
        printf("one:   K = %ld, cost = %g, cancel %g\n",                K_one, cost_one, cancellation_one);
        printf("asymp: K = %ld, cost = %g, error = %f (tol = %f)\n", K_asymp, cost_asymp, asymp_error, tolerance);
#endif

        if (asymp_error < tolerance && cost_asymp < FLINT_MIN(cost_zero, cost_one))
        {
            arb_hypgeom_legendre_p_ui_asymp(res, res_prime, n, x, K_asymp, wp);
        }
        else if (FLINT_MIN(cost_zero, cost_one) < (1e6 * prec) * prec && n < UWORD_MAX / 4)
        {
            mag_t err1, err2, xrad;
            arb_t xmid;

            mag_init(err1);
            mag_init(err2);
            mag_init(xrad);
            arb_init(xmid);

            arf_set(arb_midref(xmid), arb_midref(x));
            mag_zero(arb_radref(xmid));
            mag_set(xrad, arb_radref(x));

            arb_hypgeom_legendre_p_ui_deriv_bound(err1, err2, n, x, x2sub1);

            if (cost_zero < cost_one)
                arb_hypgeom_legendre_p_ui_zero(res, res_prime, n, xmid, K_zero, wp + cancellation_zero);
            else
                arb_hypgeom_legendre_p_ui_one(res, res_prime, n, xmid, K_one, wp + cancellation_one);

            if (res != NULL)
            {
                mag_mul(err1, err1, xrad);
                arb_add_error_mag(res, err1);
                arb_set_round(res, res, prec);
            }

            if (res_prime != NULL)
            {
                mag_mul(err2, err2, xrad);
                arb_add_error_mag(res_prime, err2);
                arb_set_round(res_prime, res_prime, prec);
            }

            mag_clear(err1);
            mag_clear(err2);
            mag_clear(xrad);
            arb_clear(xmid);
        }
        else if (asymp_error < -2.0)
        {
            /* todo -- clamp to [-1,1]? */
            arb_hypgeom_legendre_p_ui_asymp(res, res_prime, n, x, K_asymp, wp);
        }
        else
        {
            if (res != NULL)
            {
                arf_zero(arb_midref(res));
                mag_one(arb_radref(res));
            }

            if (res_prime != NULL)
            {
                arf_zero(arb_midref(res_prime));
                mag_set_ui(arb_radref(res_prime), n);
                mag_add_ui(arb_radref(res_prime), arb_radref(res_prime), 1);
                mag_mul_ui(arb_radref(res_prime), arb_radref(res_prime), n);
                mag_mul_2exp_si(arb_radref(res_prime), arb_radref(res_prime), -1);
            }
        }
    }

    arb_clear(xsub1);
    arb_clear(x2sub1);
}

