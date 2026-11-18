/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "acb.h"
#include "acb_hypgeom/sd.h"

/* The rigorous bounds of sd_bound.c against direct evaluation of the
   integrand at points of the disks, segments and rays, for parameters
   at points of the input balls. */

static double
rnd(flint_rand_t state)
{
    return n_randint(state, 20001) / 10000.0 - 1.0;
}

/* a ball, and a point of it */
static void
rand_param(acb_t x, acb_t x0, flint_rand_t state, int maxexp, int exact)
{
    double m = ldexp(1.0 + n_randint(state, 1000) / 1000.0, (int) n_randint(state, maxexp + 4) - 4);
    double th = n_randint(state, 3) ? 3.2 * rnd(state) : 1.5707963267948966 * n_randint(state, 4);

    acb_set_d_d(x, m * cos(th), n_randint(state, 4) ? m * sin(th) : 0.0);
    acb_set(x0, x);

    if (!exact && n_randint(state, 2))
    {
        mag_t r;
        double rr = m * ldexp(1.0, -1 - (int) n_randint(state, 30));

        mag_init(r);
        mag_set_d(r, rr);
        acb_add_error_mag(x, r);
        /* a point within the ball */
        arb_set_d(acb_realref(x0), arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR) + 0.7 * rr * rnd(state));
        arb_set_d(acb_imagref(x0), arf_get_d(arb_midref(acb_imagref(x)), ARF_RND_NEAR) + 0.7 * rr * rnd(state));
        mag_clear(r);
    }
}

/* Re(phi(w) - c) at the point w for the point parameters, as a double
   (accurate to an absolute error far below the tolerances used here) */
static double
re_phi(const acb_t alpha, int sym, const acb_t nu, double x, double y, const acb_t c)
{
    acb_t w, phi, d1, d2;
    double r;

    acb_init(w); acb_init(phi); acb_init(d1); acb_init(d2);
    acb_set_d_d(w, x, y);
    _sd_bound_phi(phi, d1, d2, NULL, alpha, sym, nu, w);
    acb_sub(phi, phi, c, 128);
    r = arf_get_d(arb_midref(acb_realref(phi)), ARF_RND_NEAR);
    acb_clear(w); acb_clear(phi); acb_clear(d1); acb_clear(d2);
    return r;
}

/* the weights */
static double
weight(int k, int sym, double x, double y)
{
    if (k == 0) return 1.0;
    if (k == 1) return exp(x) + (sym ? exp(-x) : 0.0);
    return hypot(x, y);
}

/* log of the mag (-inf for zero, inf for inf) */
static double
log_mag(const mag_t x)
{
    if (mag_is_zero(x)) return -D_INF;
    if (mag_is_inf(x)) return D_INF;
    return mag_get_d_log2_approx(x) * 0.6931471805599453;
}

TEST_FUNCTION_START(acb_hypgeom_sd_bound, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t alpha, alpha0, nu, nu0, c, m, a, b;
        mag_t R, M;
        mag_struct L[3];
        int sym, k, i, N;
        double mx, my, rr, ax, ay, bx, by;

        acb_init(alpha); acb_init(alpha0); acb_init(nu); acb_init(nu0);
        acb_init(c); acb_init(m); acb_init(a); acb_init(b);
        mag_init(R); mag_init(M);
        for (k = 0; k < 3; k++)
            mag_init(L + k);

        sym = n_randint(state, 2);
        rand_param(alpha, alpha0, state, 12, n_randint(state, 2));
        rand_param(nu, nu0, state, 12, n_randint(state, 2));
        if (!sym && n_randint(state, 2))
        {
            acb_set_si(alpha, -1);
            acb_set_si(alpha0, -1);
        }
        acb_set_d_d(c, 100.0 * rnd(state), 10.0 * rnd(state));

        /* the disk */
        mx = 4.0 * rnd(state);
        my = 4.0 * rnd(state);
        rr = ldexp(1.0, -(int) n_randint(state, 12)) * (1.0 + n_randint(state, 5));
        acb_set_d_d(m, mx, my);
        mag_set_d(R, rr);
        _sd_bound_disk(M, alpha, sym, nu, m, R, c);

        for (i = 0; i < 20; i++)
        {
            double th = 3.1415926535897932 * rnd(state), q = (i < 10) ? 0.999 : 0.999 * fabs(rnd(state));
            double v = re_phi(alpha0, sym, nu0, mx + q * rr * cos(th), my + q * rr * sin(th), c);

            if (v > log_mag(M) + 1e-6)
            {
                flint_printf("FAIL: disk\n\n");
                flint_printf("sym = %d, m = %g + %g i, R = %g, v = %.17g, log M = %.17g\n", sym, mx, my, rr, v, log_mag(M));
                flint_printf("alpha = "); acb_printd(alpha, 20); flint_printf("\n");
                flint_printf("nu = "); acb_printd(nu, 20); flint_printf("\n");
                flint_printf("alpha0 = "); acb_printd(alpha0, 20); flint_printf("\n");
                flint_printf("nu0 = "); acb_printd(nu0, 20); flint_printf("\n");
                flint_printf("c = "); acb_printd(c, 20); flint_printf("\n");
                flint_abort();
            }
        }

        /* the segment: the composite midpoint rule is accurate to a
           relative error of about (|phi'| h)^2 / 24 per panel, far
           below the margin, for the lengths used here */
        ax = 4.0 * rnd(state); ay = 4.0 * rnd(state);
        rr = ldexp(1.0, -(int) n_randint(state, 8));
        bx = ax + rr * rnd(state); by = ay + rr * rnd(state);
        acb_set_d_d(a, ax, ay);
        acb_set_d_d(b, bx, by);
        _sd_bound_segment(L, 3, alpha, sym, nu, a, b, c);

        N = 500;
        for (k = 0; k < 3; k++)
        {
            double sum = 0.0, top = -D_INF, len = hypot(bx - ax, by - ay);
            double * v = flint_malloc(sizeof(double) * N);

            for (i = 0; i < N; i++)
            {
                double t = (i + 0.5) / N;
                v[i] = re_phi(alpha0, sym, nu0, ax + t * (bx - ax), ay + t * (by - ay), c);
                top = FLINT_MAX(top, v[i]);
            }
            for (i = 0; i < N; i++)
            {
                double t = (i + 0.5) / N;
                sum += exp(v[i] - top) * weight(k, sym, ax + t * (bx - ax), ay + t * (by - ay));
            }
            flint_free(v);

            /* (a factor 0.98 for the quadrature) */
            if (len > 0.0 && sum > 0.0 && top + log(sum * len / N * 0.98) > log_mag(L + k) + 1e-6)
            {
                flint_printf("FAIL: segment\n\n");
                flint_printf("sym = %d, k = %d, a = %g + %g i, b = %g + %g i\n", sym, k, ax, ay, bx, by);
                flint_printf("log integral = %.17g, log bound = %.17g\n", top + log(sum * len / N), log_mag(L + k));
                flint_printf("alpha = "); acb_printd(alpha, 20); flint_printf("\n");
                flint_printf("nu = "); acb_printd(nu, 20); flint_printf("\n");
                flint_printf("alpha0 = "); acb_printd(alpha0, 20); flint_printf("\n");
                flint_printf("nu0 = "); acb_printd(nu0, 20); flint_printf("\n");
                flint_printf("c = "); acb_printd(c, 20); flint_printf("\n");
                flint_abort();
            }
            if (k == 0 && iter % 8)
                break;
        }

        /* the rays: the integral over an initial part of the ray is
           below the bound, and a horizontal ray is certified only
           inside the valley with the given index */
        {
            cd we, d;
            sd_end end;
            double ar = arf_get_d(arb_midref(acb_realref(alpha0)), ARF_RND_NEAR);
            double ai = arf_get_d(arb_midref(acb_imagref(alpha0)), ARF_RND_NEAR);
            double nr = arf_get_d(arb_midref(acb_realref(nu0)), ARF_RND_NEAR);
            double ni = arf_get_d(arb_midref(acb_imagref(nu0)), ARF_RND_NEAR);
            double tmax, off = 0.0;
            int type = (int) n_randint(state, 3), kk = (int) n_randint(state, 5) - 2;

            we.re = 6.0 * rnd(state);
            d.re = d.im = 0.0;

            if (type == 0)
            {
                end = sd_end_make(SD_PLUS, kk);
                off = 1.7 * rnd(state);
                we.im = 3.141592653589793 - atan2(ai, ar) + 6.283185307179586 * kk + off;
                d.re = 1.0;
            }
            else if (type == 1)
            {
                end = sd_end_make(SD_MINUS, kk);
                off = 1.7 * rnd(state);
                we.im = atan2(ai, ar) + 6.283185307179586 * kk + off;
                d.re = -1.0;
            }
            else
            {
                end = sd_end_make(SD_LINEAR, 0);
                we.im = 6.0 * rnd(state);
                if (sym || nr > 0.0 || n_randint(state, 2))
                    d.im = (ni > 0.0) ? -1.0 : 1.0;
                else
                {
                    d.re = nr / hypot(nr, ni);
                    d.im = -ni / hypot(nr, ni);
                }
                /* (also directions which are not valid) */
                if (n_randint(state, 10) == 0)
                    d.re = 0.3 * rnd(state);
            }

            if (_sd_bound_tail(L, 3, alpha, sym, nu, we, c, end, d))
            {
                /* far enough for the integrand to have dropped by e^-60
                   or so in the cases where the bounds are finite */
                N = 2000;
                tmax = (type == 2) ? 200.0 / (fabs(nr * d.re - ni * d.im) + 1e-3) : 12.0;
                tmax = FLINT_MIN(tmax, 300.0);

                for (k = 0; k < 3; k++)
                {
                    double sum = 0.0, top = -D_INF;
                    double * v = flint_malloc(sizeof(double) * N);
                    /* outside the valley (up to the radius of alpha) */
                    int up = (type != 2) && acb_is_exact(alpha) && fabs(off) > 1.5707963267948966 + 1e-9;

                    for (i = 0; i < N; i++)
                    {
                        double t = tmax * (i + 0.5) / N;
                        v[i] = re_phi(alpha0, sym, nu0, we.re + t * d.re, we.im + t * d.im, c);
                        top = FLINT_MAX(top, v[i]);
                    }
                    for (i = 0; i < N; i++)
                    {
                        double t = tmax * (i + 0.5) / N;
                        sum += exp(v[i] - top) * weight(k, sym, we.re + t * d.re, we.im + t * d.im);
                    }
                    flint_free(v);

                    /* (a factor 0.9 for the quadrature) */
                    if (up || (sum > 0.0 && top + log(sum * tmax / N * hypot(d.re, d.im) * 0.9) > log_mag(L + k) + 1e-6))
                    {
                        flint_printf("FAIL: tail\n\n");
                        flint_printf("sym = %d, type = %d, k = %d, kk = %d, we = %g + %g i, d = %g + %g i, up = %d\n",
                            sym, type, k, kk, we.re, we.im, d.re, d.im, up);
                        flint_printf("log integral = %.17g, log bound = %.17g\n",
                            top + log(sum * tmax / N), log_mag(L + k));
                        flint_printf("alpha = "); acb_printd(alpha, 20); flint_printf("\n");
                        flint_printf("nu = "); acb_printd(nu, 20); flint_printf("\n");
                        flint_printf("alpha0 = "); acb_printd(alpha0, 20); flint_printf("\n");
                        flint_printf("nu0 = "); acb_printd(nu0, 20); flint_printf("\n");
                        flint_printf("c = "); acb_printd(c, 20); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }

        acb_clear(alpha); acb_clear(alpha0); acb_clear(nu); acb_clear(nu0);
        acb_clear(c); acb_clear(m); acb_clear(a); acb_clear(b);
        mag_clear(R); mag_clear(M);
        for (k = 0; k < 3; k++)
            mag_clear(L + k);
    }

    TEST_FUNCTION_END(state);
}
