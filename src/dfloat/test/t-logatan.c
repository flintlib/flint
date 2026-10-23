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

/* log and atan: the static relative bounds |dN_f(x) - f(x)| <= EPS_N
   |f(x)| over the whole range (near 1 and near the table points for
   log, tiny and huge arguments, non-canonical inputs, the arb
   fallbacks); the special values; the ball versions contain f of
   sample points and are tight for exact inputs; the vector versions
   meet the same bounds. */

static void
check_rel_la(const arb_t ax, const double * y, int n, int is_log, const char * where, int which)
{
    arb_t ay, af, err, tol;
    arf_t e;
    arb_init(ay); arb_init(af); arb_init(err); arb_init(tol);
    arf_init(e);
    _dfloat_get_arb(ay, y, n, 0.0);
    if (is_log)
        arb_log(af, ax, 1500);
    else
        arb_atan(af, ax, 1500);
    arb_sub(err, ay, af, 1500);
    arb_abs(err, err);
    arb_abs(tol, af);
    arf_set_d(e, is_log ? DFLOAT_LOG_EPS(n) : DFLOAT_ATAN_EPS(n));
    arb_mul_arf(tol, tol, e, 1500);
    if (!arb_is_finite(ay) || !arb_le(err, tol))
    {
        flint_printf("FAIL: %s static bound (%s, n = %d, which = %d)\n", is_log ? "log" : "atan", where, n, which);
        flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
        flint_printf("y = "); arb_printd(ay, 40); flint_printf("\n");
        flint_printf("f = "); arb_printd(af, 40); flint_printf("\n");
        flint_printf("err = "); arb_printd(err, 10); flint_printf(", tol = "); arb_printd(tol, 10); flint_printf("\n");
        flint_abort();
    }
    arb_clear(ay); arb_clear(af); arb_clear(err); arb_clear(tol);
    arf_clear(e);
}

/* a random argument (positive for log): which 0: near 1 (1 + tiny to
   many bits, both signs); 1: near a table point 1 + i/64 or i/64; 2:
   wide range of exponents; 3: tiny; 4: huge; 5: non-canonical
   components; else moderate */
static void
random_arg_la(double * x, int n, int which, int is_log, flint_rand_t state)
{
    int k;
    if (which == 0)
    {
        _dfloat_randtest(x, n, state);
        x[0] = 1.0;
        for (k = 1; k < n; k++)
            x[k] *= d_randtest_signed(state, -300, -1);
        if (n == 1 || n_randint(state, 2))
            x[0] += d_randtest_signed(state, -60, -1);
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 1)
    {
        _dfloat_randtest(x, n, state);
        x[0] = (is_log ? 1.0 : 0.0) + (double) n_randint(state, 65) / 64.0;
        if (!is_log && n_randint(state, 2))
            x[0] = 1.0 / x[0];
        for (k = 1; k < n; k++)
            x[k] *= x[0] * d_randtest_signed(state, -300, -1);
        if (n_randint(state, 2))
            x[0] += d_randtest_signed(state, -60, -7) * x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 2)
    {
        _dfloat_randtest(x, n, state);
        x[0] = d_randtest_signed(state, -1020, 1020);
        for (k = 1; k < n; k++)
            x[k] *= x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 3)
    {
        _dfloat_randtest(x, n, state);
        x[0] = d_randtest_signed(state, -1074, -100);
        for (k = 1; k < n; k++)
            x[k] *= x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 4)
    {
        _dfloat_randtest(x, n, state);
        x[0] = d_randtest_signed(state, 100, 1023);
        for (k = 1; k < n; k++)
            x[k] *= x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 5)
    {
        for (k = 0; k < n; k++)
            x[k] = d_randtest_signed(state, -60, 9);
        if (is_log)
            x[0] = fabs(x[0]) + 1e-10;
    }
    else
    {
        _dfloat_randtest(x, n, state);
        if (n_randint(state, 2))
            x[0] = d_randtest_signed(state, -20, 10);
        _dfloat_renorm(x, n, NULL, x, n);
    }
    if (is_log && which != 5)
    {
        x[0] = fabs(x[0]);
        /* keep the value positive: the tail is below the head */
    }
    if (n_randint(state, 4) == 0 && n >= 2)
        x[n_randint(state, n)] = 0.0;
}

TEST_FUNCTION_START(logatan, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, which, is_log;
        double x[DFLOAT_MAX_N], y[DFLOAT_MAX_N];
        arb_t ax;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        which = n_randint(state, 7);
        is_log = n_randint(state, 2);
        random_arg_la(x, n, which, is_log, state);

        if (is_log)
            _dfloat_log(y, n, x);
        else
            _dfloat_atan(y, n, x);
        arb_init(ax);
        _dfloat_get_arb(ax, x, n, 0.0);

        if (arb_is_finite(ax) && (!is_log || arb_is_positive(ax)))
        {
            if (is_log && arb_is_one(ax))
            {
                if (!(y[0] == 0.0 && _dfloat_abs_sum(y, n) == 0.0))
                {
                    flint_printf("FAIL: log(1) != 0\n");
                    flint_abort();
                }
            }
            else if (!is_log && arb_is_zero(ax))
            {
                if (!(y[0] == 0.0 && _dfloat_abs_sum(y, n) == 0.0))
                {
                    flint_printf("FAIL: atan(0) != 0\n");
                    flint_abort();
                }
            }
            else
                check_rel_la(ax, y, n, is_log, "scalar", which);
        }
        else
        {
            /* special values */
            double v = y[0];
            if (is_log)
            {
                if (arb_is_zero(ax) ? !(v == -D_INF) : (arf_is_pos_inf(arb_midref(ax)) ? !(v == D_INF) : !(v != v)))
                {
                    flint_printf("FAIL: log special value (n = %d): %g\n", n, v);
                    flint_abort();
                }
            }
            else
            {
                if (arf_is_pos_inf(arb_midref(ax)) ? !(fabs(v - 1.5707963267948966) <= 0x1p-50) : (arf_is_neg_inf(arb_midref(ax)) ? !(fabs(v + 1.5707963267948966) <= 0x1p-50) : !(v != v)))
                {
                    flint_printf("FAIL: atan special value (n = %d): %g\n", n, v);
                    flint_abort();
                }
            }
        }
        arb_clear(ax);
    }

    /* balls */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, k, special, is_log, status;
        double x[DFLOAT_MAX_N + 1], y[DFLOAT_MAX_N + 1], rx, ry = 0.0;
        arb_t ax, ay, ap;
        arf_t px, t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        special = (n_randint(state, 4) == 0);
        is_log = n_randint(state, 2);

        if (special)
            _dfloat_randtest_special(x, n, state);
        else
            random_arg_la(x, n, n_randint(state, 7), is_log, state);
        rx = special ? _dfloat_randtest_rad(state) : (n_randint(state, 2) ? 0.0 : fabs(d_randtest_signed(state, -200, 2)) * (n_randint(state, 2) ? 1.0 : fabs(x[0])));

        arb_init(ax); arb_init(ay); arb_init(ap);
        arf_init(px); arf_init(t);
        _dfloat_get_arb(ax, x, n, rx);

        if (is_log)
        {
            status = _dfloat_ball_log(y, &ry, n, x, rx);
            if (status != GR_SUCCESS)
            {
                /* only allowed when the ball is not inside (0, inf) */
                if (arb_is_positive(ax) && arb_is_finite(ax))
                {
                    flint_printf("FAIL: log of a positive ball failed (n = %d)\n", n);
                    flint_printf("x = "); arb_printd(ax, 30); flint_printf(" (rad %.17g)\n", rx);
                    flint_abort();
                }
                arb_clear(ax); arb_clear(ay); arb_clear(ap);
                arf_clear(px); arf_clear(t);
                continue;
            }
        }
        else
            _dfloat_ball_atan(y, &ry, n, x, rx);
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
            if (!arf_is_finite(px) || (is_log && arf_sgn(px) <= 0))
                continue;

            arb_set_arf(ap, px);
            if (is_log)
                arb_log(ap, ap, 2000);
            else
                arb_atan(ap, ap, 2000);

            if (!arb_contains(ay, ap))
            {
                flint_printf("FAIL: %s containment\n", is_log ? "log" : "atan");
                flint_printf("n = %d, k = %d\n", n, k);
                flint_printf("x = "); arb_printd(ax, 30); flint_printf(" (rad %.17g)\n", rx);
                flint_printf("px = "); arf_printd(px, 30); flint_printf("\n");
                flint_printf("y = "); arb_printd(ay, 30); flint_printf(" (rad %.17g)\n", ry);
                flint_printf("f = "); arb_printd(ap, 30); flint_printf("\n");
                flint_abort();
            }
        }

        /* tightness for exact inputs within the kernels' ranges */
        if (rx == 0.0 && arb_is_finite(ax) && (!is_log || (x[0] >= 0x1p-1022 && x[0] < 0x1p1022)))
        {
            double ay0 = _dfloat_abs_sum(y, n);
            if (!(ry <= 2.0 * (is_log ? DFLOAT_LOG_EPS(n) : DFLOAT_ATAN_EPS(n)) * ay0 + 0x1p-1069))
            {
                flint_printf("FAIL: %s radius too large\n", is_log ? "log" : "atan");
                flint_printf("n = %d, x = %a %a %a %a, y = %a, rad = 2^%g\n", n, x[0], n > 1 ? x[1] : 0.0, n > 2 ? x[2] : 0.0, n > 3 ? x[3] : 0.0, y[0], log2(ry));
                flint_abort();
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ap);
        arf_clear(px); arf_clear(t);
    }

    /* the vector versions */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int n = 1 + n_randint(state, DFLOAT_MAX_N), k, is_log = n_randint(state, 2), status;
        slong len = n_randint(state, 12), i;
        double * x, * y, * xb, * yb;
        arb_t ax, ay, ae;
        arf_t t;

        x = flint_malloc(sizeof(double) * n * (len + 1));
        y = flint_malloc(sizeof(double) * n * (len + 1));
        xb = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
        yb = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
        arb_init(ax); arb_init(ay); arb_init(ae);
        arf_init(t);

        for (i = 0; i < len; i++)
        {
            double * xi = x + i * n, * xbi = xb + i * (n + 1);
            if (n_randint(state, 5) == 0)
                _dfloat_randtest_special(xi, n, state);
            else
                random_arg_la(xi, n, n_randint(state, 7), is_log, state);
            for (k = 0; k < n; k++)
                xbi[k] = xi[k];
            xbi[n] = n_randint(state, 3) ? 0.0 : (n_randint(state, 2) ? fabs(d_randtest_signed(state, -200, -10)) * fabs(xi[0]) : _dfloat_randtest_rad(state));
        }

        if (is_log)
        {
            switch (n)
            {
                case 1: _d1_vec_log((d1_ptr) y, (d1_srcptr) x, len); status = _d1b_vec_log((d1b_ptr) yb, (d1b_srcptr) xb, len); break;
                case 2: _d2_vec_log((d2_ptr) y, (d2_srcptr) x, len); status = _d2b_vec_log((d2b_ptr) yb, (d2b_srcptr) xb, len); break;
                case 3: _d3_vec_log((d3_ptr) y, (d3_srcptr) x, len); status = _d3b_vec_log((d3b_ptr) yb, (d3b_srcptr) xb, len); break;
                default: _d4_vec_log((d4_ptr) y, (d4_srcptr) x, len); status = _d4b_vec_log((d4b_ptr) yb, (d4b_srcptr) xb, len); break;
            }
        }
        else
        {
            status = GR_SUCCESS;
            switch (n)
            {
                case 1: _d1_vec_atan((d1_ptr) y, (d1_srcptr) x, len); _d1b_vec_atan((d1b_ptr) yb, (d1b_srcptr) xb, len); break;
                case 2: _d2_vec_atan((d2_ptr) y, (d2_srcptr) x, len); _d2b_vec_atan((d2b_ptr) yb, (d2b_srcptr) xb, len); break;
                case 3: _d3_vec_atan((d3_ptr) y, (d3_srcptr) x, len); _d3b_vec_atan((d3b_ptr) yb, (d3b_srcptr) xb, len); break;
                default: _d4_vec_atan((d4_ptr) y, (d4_srcptr) x, len); _d4b_vec_atan((d4b_ptr) yb, (d4b_srcptr) xb, len); break;
            }
        }

        for (i = 0; i < len; i++)
        {
            double * xi = x + i * n, * yi = y + i * n, * xbi = xb + i * (n + 1), * ybi = yb + i * (n + 1);
            double z[DFLOAT_MAX_N];

            _dfloat_get_arb(ax, xi, n, 0.0);
            if (arb_is_finite(ax) && (!is_log || arb_is_positive(ax)) && !(is_log && arb_is_one(ax)) && !(!is_log && arb_is_zero(ax)))
                check_rel_la(ax, yi, n, is_log, "vector", 0);
            else
            {
                /* the same as the scalar function */
                if (is_log)
                    _dfloat_log(z, n, xi);
                else
                    _dfloat_atan(z, n, xi);
                for (k = 0; k < n; k++)
                    if (z[k] != yi[k] && !(z[k] != z[k] && yi[k] != yi[k]))
                    {
                        flint_printf("FAIL: vector %s special value (n = %d)\n", is_log ? "log" : "atan", n);
                        flint_abort();
                    }
            }

            /* ball: containment at the midpoint and endpoints (a log
               ball not inside (0, inf) has failed as a whole) */
            _dfloat_get_arb(ax, xbi, n, xbi[n]);
            if (is_log && (status != GR_SUCCESS))
                continue;
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
                if (!arf_is_finite(t) || (is_log && arf_sgn(t) <= 0))
                    continue;
                arb_set_arf(ae, t);
                if (is_log)
                    arb_log(ae, ae, 2000);
                else
                    arb_atan(ae, ae, 2000);
                if (!arb_contains(ay, ae))
                {
                    flint_printf("FAIL: vector ball containment (%s, n = %d, i = %wd, k = %d)\n", is_log ? "log" : "atan", n, i, k);
                    flint_printf("x = "); arb_printd(ax, 30); flint_printf("\n");
                    flint_printf("y = "); arb_printd(ay, 30); flint_printf("\n");
                    flint_printf("f = "); arb_printd(ae, 30); flint_printf("\n");
                    flint_abort();
                }
            }
        }
        if (is_log && status != GR_SUCCESS)
        {
            /* some element must not be inside (0, inf) */
            int found = 0;
            for (i = 0; i < len; i++)
            {
                _dfloat_get_arb(ax, xb + i * (n + 1), n, xb[i * (n + 1) + n]);
                if (!(arb_is_positive(ax) && arb_is_finite(ax)))
                    found = 1;
            }
            if (!found)
            {
                flint_printf("FAIL: vector ball log failed on positive balls\n");
                flint_abort();
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ae);
        arf_clear(t);
        flint_free(x); flint_free(y); flint_free(xb); flint_free(yb);
    }

    TEST_FUNCTION_END(state);
}
