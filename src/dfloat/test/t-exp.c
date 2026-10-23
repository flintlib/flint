/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "dfloat.h"

/* exp on balls: containment of exp of sample points (arb at high
   precision), over the whole range including overflow and underflow;
   the radius must be reasonably tight for exact inputs of moderate
   size; exp(0) is exactly 1. */

TEST_FUNCTION_START(exp, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    /* the static bound of the plain kernel: |dN_exp(x) - exp(x)| <=
       EPS_N |exp(x)| (+ 2^-1070) over the whole domain, including
       arguments near multiples of log 2 (cancellation in the
       reduction), near the ends of the range, tiny, and with
       non-canonical components */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, k, which;
        double x[DFLOAT_MAX_N], y[DFLOAT_MAX_N];
        arb_t ax, ay, ae, err, tol;
        arf_t t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        which = n_randint(state, 6);
        arb_init(ax); arb_init(ay); arb_init(ae); arb_init(err); arb_init(tol);
        arf_init(t);

        if (which == 0)
        {
            /* x = k log 2 + tiny, to many bits */
            slong kk = (slong) n_randint(state, 2000) - 1000;
            arb_const_log2(ax, 400);
            arb_mul_si(ax, ax, kk, 400);
            arf_set_d(t, d_randtest_signed(state, -300, -20));
            arb_add_arf(ax, ax, t, 400);
            _dfloat_set_arf(x, n, NULL, arb_midref(ax));
        }
        else if (which == 1)
        {
            /* near the ends of the range */
            x[0] = n_randint(state, 2) ? 709.0 + n_randint(state, 1000) * 1e-3 : -745.0 - n_randint(state, 1000) * 1e-3;
            for (k = 1; k < n; k++)
                x[k] = 0.0;
            _dfloat_randtest(y, n, state);
            if (fabs(y[0]) < 1.0)
            {
                double e = 0.0;
                _dfloat_add_n(x, &e, n, x, y);
            }
        }
        else if (which == 2)
        {
            /* tiny */
            _dfloat_randtest(x, n, state);
            x[0] = d_randtest_signed(state, -1000, -15);
            for (k = 1; k < n; k++)
                x[k] *= x[0];
            _dfloat_renorm(x, n, NULL, x, n);
        }
        else if (which == 3)
        {
            /* non-canonical components */
            for (k = 0; k < n; k++)
                x[k] = d_randtest_signed(state, -60, 9);
        }
        else
        {
            _dfloat_randtest(x, n, state);
            if (n_randint(state, 2))
                x[0] = d_randtest_signed(state, -20, 10);
            _dfloat_renorm(x, n, NULL, x, n);
        }
        if (n_randint(state, 4) == 0 && n >= 2)
            x[n_randint(state, n)] = 0.0;

        _dfloat_exp(y, n, x);
        _dfloat_get_arb(ax, x, n, 0.0);
        _dfloat_get_arb(ay, y, n, 0.0);

        if (arb_is_finite(ax) && arf_cmp_d(arb_midref(ax), 709.9) <= 0 && arf_cmp_d(arb_midref(ax), -746.0) >= 0
            && arb_is_finite(ay))
        {
            arb_exp(ae, ax, 1500);
            arb_sub(err, ay, ae, 1500);
            arb_abs(err, err);
            arb_abs(tol, ae);
            { arf_t ee; arf_init(ee); arf_set_d(ee, DFLOAT_EXP_EPS(n)); arb_mul_arf(tol, tol, ee, 1500); arf_clear(ee); }
            arf_set_d(t, 0x1p-1070);
            arb_add_arf(tol, tol, t, 1500);
            if (!arb_lt(err, tol))
            {
                flint_printf("FAIL: static bound (n = %d, which = %d)\n", n, which);
                flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
                flint_printf("y = "); arb_printd(ay, 40); flint_printf("\n");
                flint_printf("exp = "); arb_printd(ae, 40); flint_printf("\n");
                flint_printf("err = "); arb_printd(err, 10); flint_printf(", tol = "); arb_printd(tol, 10); flint_printf("\n");
                flint_abort();
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ae); arb_clear(err); arb_clear(tol);
        arf_clear(t);
    }

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, k, special;
        double x[DFLOAT_MAX_N + 1], y[DFLOAT_MAX_N + 1], rx, ry;
        arb_t ax, ay, ap;
        arf_t px, t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        special = (n_randint(state, 4) == 0);

        if (special)
            _dfloat_randtest_special(x, n, state);
        else
        {
            _dfloat_randtest(x, n, state);
            /* moderate arguments most of the time */
            if (n_randint(state, 2))
                x[0] = d_randtest_signed(state, -20, 10);
            if (n_randint(state, 2))
                for (k = 1; k < n; k++)
                    x[k] = 0.0;
            _dfloat_renorm(x, n, NULL, x, n);
        }
        rx = special ? _dfloat_randtest_rad(state) : (n_randint(state, 2) ? 0.0 : fabs(d_randtest_signed(state, -200, -20)));

        _dfloat_ball_exp(y, &ry, n, x, rx);

        arb_init(ax); arb_init(ay); arb_init(ap);
        arf_init(px); arf_init(t);
        _dfloat_get_arb(ax, x, n, rx);
        _dfloat_get_arb(ay, y, n, ry);

        for (k = 0; k < 4; k++)
        {
            /* sample points: midpoint, endpoints, random */
            if (k == 0 || rx == D_INF)
            {
                arf_set(px, arb_midref(ax));
                if (k > 0)
                {
                    arf_set_d(t, d_randtest_signed(state, -1074, 1023));
                    arf_add(px, px, t, ARF_PREC_EXACT, ARF_RND_DOWN);
                }
            }
            else
            {
                arf_set_d(t, rx);
                if (k == 3)
                {
                    arf_mul_ui(t, t, n_randint(state, 1000), 30, ARF_RND_DOWN);
                    arf_mul_2exp_si(t, t, -10);
                }
                if (k == 2 || (k == 3 && n_randint(state, 2)))
                    arf_neg(t, t);
                arf_add(px, arb_midref(ax), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            if (!arf_is_finite(px))
                continue;

            arb_set_arf(ap, px);
            arb_exp(ap, ap, 2000);

            if (!arb_contains(ay, ap))
            {
                flint_printf("FAIL: containment\n");
                flint_printf("n = %d, k = %d\n", n, k);
                flint_printf("x = "); arb_printd(ax, 30); flint_printf(" (rad %.17g)\n", rx);
                flint_printf("px = "); arf_printd(px, 30); flint_printf("\n");
                flint_printf("y = "); arb_printd(ay, 30); flint_printf(" (rad %.17g)\n", ry);
                flint_printf("exp = "); arb_printd(ap, 30); flint_printf("\n");
                flint_abort();
            }
        }

        /* tightness for exact moderate inputs */
        if (rx == 0.0 && fabs(x[0]) < 500.0 && ry != D_INF)
        {
            double rel = ry / fabs(y[0]);
            if (!(rel < 4.0 * DFLOAT_EXP_EPS(n)))
            {
                flint_printf("FAIL: radius too large\n");
                flint_printf("n = %d, x = %a %a %a %a, rel rad = 2^%g\n", n, x[0], n > 1 ? x[1] : 0.0, n > 2 ? x[2] : 0.0, n > 3 ? x[3] : 0.0, log2(rel));
                flint_abort();
            }
        }

        if (x[0] == 0.0 && rx == 0.0 && _dfloat_abs_sum(x, n) == 0.0)
        {
            if (!(y[0] == 1.0 && ry == 0.0))
            {
                flint_printf("FAIL: exp(0) != 1 exactly\n");
                flint_abort();
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ap);
        arf_clear(px); arf_clear(t);
    }

    /* the vector versions: the plain one meets the static bound on
       every lane, the ball one contains the scalar ball's sample
       points (exp of the midpoint and of the endpoints) */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int n = 1 + n_randint(state, DFLOAT_MAX_N), k;
        slong len = n_randint(state, 12), i;
        double * x, * y, * xb, * yb;
        arb_t ax, ay, ae, err, tol;
        arf_t t;

        x = flint_malloc(sizeof(double) * n * (len + 1));
        y = flint_malloc(sizeof(double) * n * (len + 1));
        xb = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
        yb = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
        arb_init(ax); arb_init(ay); arb_init(ae); arb_init(err); arb_init(tol);
        arf_init(t);

        for (i = 0; i < len; i++)
        {
            double * xi = x + i * n, * xbi = xb + i * (n + 1);
            switch (n_randint(state, 4))
            {
                case 0: _dfloat_randtest_special(xi, n, state); break;
                case 1:
                    _dfloat_randtest(xi, n, state);
                    xi[0] = d_randtest_signed(state, -20, 10);
                    _dfloat_renorm(xi, n, NULL, xi, n);
                    break;
                case 2:
                    for (k = 0; k < n; k++)
                        xi[k] = d_randtest_signed(state, -60, 9);
                    break;
                default:
                    _dfloat_randtest(xi, n, state);
                    if (n_randint(state, 2))
                        xi[0] = (double) ((slong) n_randint(state, 1500) - 750);
                    _dfloat_renorm(xi, n, NULL, xi, n);
            }
            for (k = 0; k < n; k++)
                xbi[k] = xi[k];
            xbi[n] = n_randint(state, 3) ? 0.0 : (n_randint(state, 2) ? fabs(d_randtest_signed(state, -200, -20)) : _dfloat_randtest_rad(state));
        }

        switch (n)
        {
            case 1: _d1_vec_exp((d1_ptr) y, (d1_srcptr) x, len); _d1b_vec_exp((d1b_ptr) yb, (d1b_srcptr) xb, len); break;
            case 2: _d2_vec_exp((d2_ptr) y, (d2_srcptr) x, len); _d2b_vec_exp((d2b_ptr) yb, (d2b_srcptr) xb, len); break;
            case 3: _d3_vec_exp((d3_ptr) y, (d3_srcptr) x, len); _d3b_vec_exp((d3b_ptr) yb, (d3b_srcptr) xb, len); break;
            default: _d4_vec_exp((d4_ptr) y, (d4_srcptr) x, len); _d4b_vec_exp((d4b_ptr) yb, (d4b_srcptr) xb, len); break;
        }

        for (i = 0; i < len; i++)
        {
            double * xi = x + i * n, * yi = y + i * n, * xbi = xb + i * (n + 1), * ybi = yb + i * (n + 1);
            double z[DFLOAT_MAX_N];

            /* plain: the static bound, or the same as the scalar
               function outside the domain */
            _dfloat_get_arb(ax, xi, n, 0.0);
            _dfloat_get_arb(ay, yi, n, 0.0);
            if (arb_is_finite(ax) && arf_cmp_d(arb_midref(ax), 709.9) <= 0 && arf_cmp_d(arb_midref(ax), -746.0) >= 0)
            {
                arb_exp(ae, ax, 1500);
                arb_sub(err, ay, ae, 1500);
                arb_abs(err, err);
                arb_abs(tol, ae);
                { arf_t ee; arf_init(ee); arf_set_d(ee, DFLOAT_EXP_EPS(n)); arb_mul_arf(tol, tol, ee, 1500); arf_clear(ee); }
                arf_set_d(t, 0x1p-1070);
                arb_add_arf(tol, tol, t, 1500);
                if (!arb_is_finite(ay) || !arb_lt(err, tol))
                {
                    flint_printf("FAIL: vector static bound (n = %d, i = %wd)\n", n, i);
                    flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
                    flint_printf("y = "); arb_printd(ay, 40); flint_printf("\n");
                    flint_printf("exp = "); arb_printd(ae, 40); flint_printf("\n");
                    flint_abort();
                }
            }
            else
            {
                _dfloat_exp(z, n, xi);
                for (k = 0; k < n; k++)
                    if (z[k] != yi[k] && !(z[k] != z[k] && yi[k] != yi[k]))
                    {
                        flint_printf("FAIL: vector exp outside the domain (n = %d)\n", n);
                        flint_abort();
                    }
            }

            /* ball: containment of exp at the midpoint and endpoints */
            _dfloat_get_arb(ax, xbi, n, xbi[n]);
            _dfloat_get_arb(ay, ybi, n, ybi[n]);
            for (k = 0; k < 3; k++)
            {
                arf_set(t, arb_midref(ax));
                if (k > 0 && xbi[n] != D_INF)
                {
                    arf_t rr;
                    arf_init(rr);
                    arf_set_d(rr, k == 1 ? xbi[n] : -xbi[n]);
                    arf_add(t, t, rr, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_clear(rr);
                }
                if (!arf_is_finite(t))
                    continue;
                arb_set_arf(ae, t);
                arb_exp(ae, ae, 2000);
                if (!arb_contains(ay, ae))
                {
                    flint_printf("FAIL: vector ball containment (n = %d, i = %wd, k = %d)\n", n, i, k);
                    flint_printf("x = "); arb_printd(ax, 30); flint_printf("\n");
                    flint_printf("y = "); arb_printd(ay, 30); flint_printf("\n");
                    flint_printf("exp = "); arb_printd(ae, 30); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ae); arb_clear(err); arb_clear(tol);
        arf_clear(t);
        flint_free(x); flint_free(y); flint_free(xb); flint_free(yb);
    }

    TEST_FUNCTION_END(state);
}
