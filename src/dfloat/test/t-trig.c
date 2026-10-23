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

/* sin, cos, sin_cos: the static bounds |dN_sin(x) - sin x| <= EPS_N
   (|x| + min(|x|, 1)) and |dN_cos(x) - cos x| <= EPS_N (|x| + 1) over
   the whole range (near multiples of pi/2, large, tiny, non-canonical
   inputs, and beyond the kernel's range where arb is used); sin_cos
   agrees with sin and cos; the ball versions contain sin and cos of
   sample points and are tight for exact inputs; the vector versions
   meet the same bounds. */

/* the tolerance EPS_n (|x| + a) as an arb */
static void
trig_tol(arb_t tol, const arb_t ax, int n, int is_sin)
{
    arb_t a;
    arf_t e;
    arb_init(a);
    arf_init(e);
    arb_abs(tol, ax);
    if (is_sin)
    {
        arb_one(a);
        arb_min(a, a, tol, 1500);
    }
    else
        arb_one(a);
    arb_add(tol, tol, a, 1500);
    arf_set_d(e, DFLOAT_TRIG_EPS(n));
    arb_mul_arf(tol, tol, e, 1500);
    arb_clear(a);
    arf_clear(e);
}

static void
check_static(const arb_t ax, const double * y, int n, int is_sin, const char * where, int which)
{
    arb_t ay, ae, err, tol;
    arb_init(ay); arb_init(ae); arb_init(err); arb_init(tol);
    _dfloat_get_arb(ay, y, n, 0.0);
    if (is_sin)
        arb_sin(ae, ax, 1500);
    else
        arb_cos(ae, ax, 1500);
    arb_sub(err, ay, ae, 1500);
    arb_abs(err, err);
    trig_tol(tol, ax, n, is_sin);
    if (!arb_is_finite(ay) || !arb_le(err, tol))
    {
        flint_printf("FAIL: %s static bound (%s, n = %d, which = %d)\n", is_sin ? "sin" : "cos", where, n, which);
        flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
        flint_printf("y = "); arb_printd(ay, 40); flint_printf("\n");
        flint_printf("f = "); arb_printd(ae, 40); flint_printf("\n");
        flint_printf("err = "); arb_printd(err, 10); flint_printf(", tol = "); arb_printd(tol, 10); flint_printf("\n");
        flint_abort();
    }
    arb_clear(ay); arb_clear(ae); arb_clear(err); arb_clear(tol);
}

/* a random argument: which 0: k pi/2 + tiny to many bits (|k| up to
   2^22, beyond the kernel's range); 1: large, up to and beyond the
   range; 2: tiny; 3: non-canonical components; else moderate */
static void
random_arg(double * x, int n, int which, flint_rand_t state)
{
    int k;
    if (which == 0)
    {
        arb_t ax;
        arf_t t;
        slong kk = (slong) n_randint(state, 1 << 23) - (1 << 22);
        arb_init(ax);
        arf_init(t);
        arb_const_pi(ax, 400);
        arb_mul_2exp_si(ax, ax, -1);
        arb_mul_si(ax, ax, kk, 400);
        if (n_randint(state, 4))
        {
            arf_set_d(t, d_randtest_signed(state, -300, -20));
            arb_add_arf(ax, ax, t, 400);
        }
        _dfloat_set_arf(x, n, NULL, arb_midref(ax));
        arb_clear(ax);
        arf_clear(t);
    }
    else if (which == 1)
    {
        _dfloat_randtest(x, n, state);
        x[0] = d_randtest_signed(state, 15, 30);
        for (k = 1; k < n; k++)
            x[k] *= x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 2)
    {
        _dfloat_randtest(x, n, state);
        x[0] = d_randtest_signed(state, -1000, -10);
        for (k = 1; k < n; k++)
            x[k] *= x[0];
        _dfloat_renorm(x, n, NULL, x, n);
    }
    else if (which == 3)
    {
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
}

TEST_FUNCTION_START(trig, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, k, which, want;
        double x[DFLOAT_MAX_N], s[DFLOAT_MAX_N], c[DFLOAT_MAX_N], s2[DFLOAT_MAX_N], c2[DFLOAT_MAX_N];
        arb_t ax;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        which = n_randint(state, 6);
        want = 1 + n_randint(state, 3);
        random_arg(x, n, which, state);

        _dfloat_sin_cos(s, c, n, x, want);
        arb_init(ax);
        _dfloat_get_arb(ax, x, n, 0.0);

        if (arb_is_finite(ax))
        {
            if (want & 1)
                check_static(ax, s, n, 1, "scalar", which);
            if (want & 2)
                check_static(ax, c, n, 0, "scalar", which);

            /* sin_cos agrees with sin and cos */
            if (want == 3)
            {
                _dfloat_sin_cos(s2, c2, n, x, 1);
                _dfloat_sin_cos(s2, c2, n, x, 2);
                for (k = 0; k < n; k++)
                {
                    if (s[k] != s2[k] || c[k] != c2[k])
                    {
                        flint_printf("FAIL: sin_cos vs sin, cos (n = %d, which = %d)\n", n, which);
                        flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
                        flint_printf("%a %a\n", s[k], s2[k]);
                        flint_printf("%a %a\n", c[k], c2[k]);
                        flint_abort();
                    }
                }
            }

            /* exact zero */
            if (_dfloat_abs_sum(x, n) == 0.0)
            {
                if (((want & 1) && (s[0] != 0.0 || _dfloat_abs_sum(s, n) != 0.0))
                    || ((want & 2) && (c[0] != 1.0 || _dfloat_abs_sum(c + 1, n - 1) != 0.0)))
                {
                    flint_printf("FAIL: sin(0), cos(0)\n");
                    flint_abort();
                }
            }
        }
        else
        {
            if (((want & 1) && !(s[0] != s[0])) || ((want & 2) && !(c[0] != c[0])))
            {
                flint_printf("FAIL: non-finite input should give nan\n");
                flint_abort();
            }
        }
        arb_clear(ax);
    }

    /* balls */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, k, special, want, f;
        double x[DFLOAT_MAX_N + 1], s[DFLOAT_MAX_N + 1], c[DFLOAT_MAX_N + 1], rx, rs, rc;
        arb_t ax, ay, ap;
        arf_t px, t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        special = (n_randint(state, 4) == 0);
        want = 1 + n_randint(state, 3);

        if (special)
            _dfloat_randtest_special(x, n, state);
        else
            random_arg(x, n, n_randint(state, 6), state);
        rx = special ? _dfloat_randtest_rad(state) : (n_randint(state, 2) ? 0.0 : fabs(d_randtest_signed(state, -200, 2)));

        _dfloat_ball_sin_cos(s, &rs, c, &rc, n, x, rx, want);

        arb_init(ax); arb_init(ay); arb_init(ap);
        arf_init(px); arf_init(t);
        _dfloat_get_arb(ax, x, n, rx);

        for (f = 0; f < 2; f++)
        {
            if (!(want & (1 << f)))
                continue;
            _dfloat_get_arb(ay, f == 0 ? s : c, n, f == 0 ? rs : rc);

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
                if (f == 0)
                    arb_sin(ap, ap, 2000);
                else
                    arb_cos(ap, ap, 2000);

                if (!arb_contains(ay, ap))
                {
                    flint_printf("FAIL: %s containment\n", f == 0 ? "sin" : "cos");
                    flint_printf("n = %d, k = %d, want = %d\n", n, k, want);
                    flint_printf("x = "); arb_printd(ax, 30); flint_printf(" (rad %.17g)\n", rx);
                    flint_printf("px = "); arf_printd(px, 30); flint_printf("\n");
                    flint_printf("y = "); arb_printd(ay, 30); flint_printf(" (rad %.17g)\n", f == 0 ? rs : rc);
                    flint_printf("f = "); arb_printd(ap, 30); flint_printf("\n");
                    flint_abort();
                }
            }

            /* tightness for exact inputs within the kernel's range */
            if (rx == 0.0 && arb_is_finite(ax) && fabs(x[0]) <= 1.6e6)
            {
                double ry = (f == 0) ? rs : rc;
                double am = _dfloat_abs_sum(x, n);
                if (!(ry <= 2.0 * DFLOAT_TRIG_EPS(n) * (am + 1.0) + 0x1p-1069))
                {
                    flint_printf("FAIL: %s radius too large\n", f == 0 ? "sin" : "cos");
                    flint_printf("n = %d, x = %a %a %a %a, rad = 2^%g\n", n, x[0], n > 1 ? x[1] : 0.0, n > 2 ? x[2] : 0.0, n > 3 ? x[3] : 0.0, log2(ry));
                    flint_abort();
                }
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ap);
        arf_clear(px); arf_clear(t);
    }

    /* the vector versions: the plain ones meet the static bounds on
       every lane, the ball ones contain sin / cos of the midpoint and
       endpoints */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int n = 1 + n_randint(state, DFLOAT_MAX_N), k, f;
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
                random_arg(xi, n, n_randint(state, 6), state);
            for (k = 0; k < n; k++)
                xbi[k] = xi[k];
            xbi[n] = n_randint(state, 3) ? 0.0 : (n_randint(state, 2) ? fabs(d_randtest_signed(state, -200, 2)) : _dfloat_randtest_rad(state));
        }

        for (f = 0; f < 2; f++)
        {
            if (f == 0)
            {
                switch (n)
                {
                    case 1: _d1_vec_sin((d1_ptr) y, (d1_srcptr) x, len); _d1b_vec_sin((d1b_ptr) yb, (d1b_srcptr) xb, len); break;
                    case 2: _d2_vec_sin((d2_ptr) y, (d2_srcptr) x, len); _d2b_vec_sin((d2b_ptr) yb, (d2b_srcptr) xb, len); break;
                    case 3: _d3_vec_sin((d3_ptr) y, (d3_srcptr) x, len); _d3b_vec_sin((d3b_ptr) yb, (d3b_srcptr) xb, len); break;
                    default: _d4_vec_sin((d4_ptr) y, (d4_srcptr) x, len); _d4b_vec_sin((d4b_ptr) yb, (d4b_srcptr) xb, len); break;
                }
            }
            else
            {
                switch (n)
                {
                    case 1: _d1_vec_cos((d1_ptr) y, (d1_srcptr) x, len); _d1b_vec_cos((d1b_ptr) yb, (d1b_srcptr) xb, len); break;
                    case 2: _d2_vec_cos((d2_ptr) y, (d2_srcptr) x, len); _d2b_vec_cos((d2b_ptr) yb, (d2b_srcptr) xb, len); break;
                    case 3: _d3_vec_cos((d3_ptr) y, (d3_srcptr) x, len); _d3b_vec_cos((d3b_ptr) yb, (d3b_srcptr) xb, len); break;
                    default: _d4_vec_cos((d4_ptr) y, (d4_srcptr) x, len); _d4b_vec_cos((d4b_ptr) yb, (d4b_srcptr) xb, len); break;
                }
            }

            for (i = 0; i < len; i++)
            {
                double * xi = x + i * n, * yi = y + i * n, * xbi = xb + i * (n + 1), * ybi = yb + i * (n + 1);

                _dfloat_get_arb(ax, xi, n, 0.0);
                if (arb_is_finite(ax))
                    check_static(ax, yi, n, f == 0, "vector", 0);
                else if (!(yi[0] != yi[0]))
                {
                    flint_printf("FAIL: vector non-finite input should give nan\n");
                    flint_abort();
                }

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
                    if (f == 0)
                        arb_sin(ae, ae, 2000);
                    else
                        arb_cos(ae, ae, 2000);
                    if (!arb_contains(ay, ae))
                    {
                        flint_printf("FAIL: vector ball containment (%s, n = %d, i = %wd, k = %d)\n", f == 0 ? "sin" : "cos", n, i, k);
                        flint_printf("x = "); arb_printd(ax, 30); flint_printf("\n");
                        flint_printf("y = "); arb_printd(ay, 30); flint_printf("\n");
                        flint_printf("f = "); arb_printd(ae, 30); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }

        /* vec_sin_cos: identical to vec_sin and vec_cos */
        {
            double * s1 = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
            double * c1 = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
            double * s2 = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
            double * c2 = flint_malloc(sizeof(double) * (n + 1) * (len + 1));
            int b;
            for (b = 0; b <= 1; b++)
            {
                slong w = n + b;
#define VSC(X) \
                if (b) { _##X##b_vec_sin((X##b_ptr) s1, (X##b_srcptr) xb, len); _##X##b_vec_cos((X##b_ptr) c1, (X##b_srcptr) xb, len); \
                         _##X##b_vec_sin_cos((X##b_ptr) s2, (X##b_ptr) c2, (X##b_srcptr) xb, len); } \
                else { _##X##_vec_sin((X##_ptr) s1, (X##_srcptr) x, len); _##X##_vec_cos((X##_ptr) c1, (X##_srcptr) x, len); \
                       _##X##_vec_sin_cos((X##_ptr) s2, (X##_ptr) c2, (X##_srcptr) x, len); }
                switch (n)
                {
                    case 1: VSC(d1) break;
                    case 2: VSC(d2) break;
                    case 3: VSC(d3) break;
                    default: VSC(d4) break;
                }
#undef VSC
                for (i = 0; i < len * w; i++)
                {
                    if ((s1[i] != s2[i] && !(s1[i] != s1[i] && s2[i] != s2[i]))
                        || (c1[i] != c2[i] && !(c1[i] != c1[i] && c2[i] != c2[i])))
                    {
                        flint_printf("FAIL: vec_sin_cos (n = %d, ball = %d, i = %wd)\n", n, b, i);
                        flint_abort();
                    }
                }
            }
            flint_free(s1); flint_free(c1); flint_free(s2); flint_free(c2);
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ae);
        arf_clear(t);
        flint_free(x); flint_free(y); flint_free(xb); flint_free(yb);
    }

    TEST_FUNCTION_END(state);
}
