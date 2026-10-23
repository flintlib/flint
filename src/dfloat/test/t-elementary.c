/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include <string.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "dfloat.h"

/* expm1 and log1p: the static relative bounds; the composed
   functions tan, sinh, cosh, tanh, asin, acos, atan2, pow: a loose
   relative bound (2^(14 - 53 n)) against arb, and the ball versions
   contain the true values at sample points. */

#define NFUN 10
static const char * fnames[NFUN] = {"expm1", "log1p", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan2", "pow"};

/* the plain function; returns 0 if undefined for the arguments (the
   caller skips the check) */
static int
eval_plain(double * res, int f, int n, const double * x, const double * y)
{
    switch (n)
    {
#define CASES(X) \
        switch (f) { \
            case 0: X##_expm1((X##_ptr) res, (X##_srcptr) x); break; \
            case 1: X##_log1p((X##_ptr) res, (X##_srcptr) x); break; \
            case 2: X##_tan((X##_ptr) res, (X##_srcptr) x); break; \
            case 3: X##_sinh((X##_ptr) res, (X##_srcptr) x); break; \
            case 4: X##_cosh((X##_ptr) res, (X##_srcptr) x); break; \
            case 5: X##_tanh((X##_ptr) res, (X##_srcptr) x); break; \
            case 6: X##_asin((X##_ptr) res, (X##_srcptr) x); break; \
            case 7: X##_acos((X##_ptr) res, (X##_srcptr) x); break; \
            case 8: X##_atan2((X##_ptr) res, (X##_srcptr) y, (X##_srcptr) x); break; \
            default: X##_pow((X##_ptr) res, (X##_srcptr) x, (X##_srcptr) y); break; }
        case 1: CASES(d1) break;
        case 2: CASES(d2) break;
        case 3: CASES(d3) break;
        default: CASES(d4) break;
#undef CASES
    }
    return 1;
}

static int
eval_ball(double * res, int f, int n, const double * x, const double * y)
{
    int status = GR_SUCCESS;
    switch (n)
    {
#define CASES(X) \
        switch (f) { \
            case 0: X##_expm1((X##_ptr) res, (X##_srcptr) x); break; \
            case 1: status = X##_log1p((X##_ptr) res, (X##_srcptr) x); break; \
            case 2: status = X##_tan((X##_ptr) res, (X##_srcptr) x); break; \
            case 3: X##_sinh((X##_ptr) res, (X##_srcptr) x); break; \
            case 4: X##_cosh((X##_ptr) res, (X##_srcptr) x); break; \
            case 5: X##_tanh((X##_ptr) res, (X##_srcptr) x); break; \
            case 6: status = X##_asin((X##_ptr) res, (X##_srcptr) x); break; \
            case 7: status = X##_acos((X##_ptr) res, (X##_srcptr) x); break; \
            case 8: X##_atan2((X##_ptr) res, (X##_srcptr) y, (X##_srcptr) x); break; \
            default: status = X##_pow((X##_ptr) res, (X##_srcptr) x, (X##_srcptr) y); break; }
        case 1: CASES(d1b) break;
        case 2: CASES(d2b) break;
        case 3: CASES(d3b) break;
        default: CASES(d4b) break;
#undef CASES
    }
    return status;
}

/* the reference; returns 0 if arb cannot evaluate (outside the domain) */
static int
eval_arb(arb_t res, int f, const arb_t x, const arb_t y, slong prec)
{
    switch (f)
    {
        case 0: arb_expm1(res, x, prec); break;
        case 1: arb_log1p(res, x, prec); break;
        case 2: arb_tan(res, x, prec); break;
        case 3: arb_sinh(res, x, prec); break;
        case 4: arb_cosh(res, x, prec); break;
        case 5: arb_tanh(res, x, prec); break;
        case 6: arb_asin(res, x, prec); break;
        case 7: arb_acos(res, x, prec); break;
        case 8: arb_atan2(res, y, x, prec); break;
        default: arb_pow(res, x, y, prec); break;
    }
    return arb_is_finite(res);
}

static void
random_arg_el(double * x, int n, int f, flint_rand_t state)
{
    int k;
    _dfloat_randtest(x, n, state);
    switch (n_randint(state, 4))
    {
        case 0: x[0] = d_randtest_signed(state, -300, -1); break;
        case 1: x[0] = ((double) n_randint(state, 2000001) - 1000000) * 1e-6; break;   /* [-1, 1] */
        case 2: x[0] = ((double) n_randint(state, 2000001) - 1000000) * 1e-4; break;   /* [-100, 100] */
        default: x[0] = d_randtest_signed(state, -20, 10); break;
    }
    for (k = 1; k < n; k++)
        x[k] *= x[0];
    _dfloat_renorm(x, n, NULL, x, n);
    if (f == 6 || f == 7)
    {
        /* asin, acos: mostly in [-1, 1], sometimes exactly +-1 */
        if (n_randint(state, 8) == 0)
        {
            x[0] = n_randint(state, 2) ? 1.0 : -1.0;
            for (k = 1; k < n; k++)
                x[k] = 0.0;
        }
        else if (fabs(x[0]) > 1.0)
        {
            double s = 1.0 / (2.0 * fabs(x[0]));
            for (k = 0; k < n; k++)
                x[k] *= s;
            _dfloat_renorm(x, n, NULL, x, n);
        }
    }
    if (f == 9)
        x[0] = fabs(x[0]);    /* pow: a positive base */
    if (n_randint(state, 4) == 0 && n >= 2)
        x[n_randint(state, n)] = 0.0;
}

TEST_FUNCTION_START(elementary, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 30000 * flint_test_multiplier(); iter++)
    {
        int n, f, k;
        double x[DFLOAT_MAX_N], y[DFLOAT_MAX_N], z[DFLOAT_MAX_N];
        arb_t ax, ay, az, af, err, tol;
        arf_t t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        f = n_randint(state, NFUN);
        random_arg_el(x, n, f, state);
        random_arg_el(y, n, f == 9 ? 0 : f, state);
        if (f == 9 && n_randint(state, 3) == 0)
        {
            /* an integer or half-integer exponent */
            y[0] = (double) ((slong) n_randint(state, 41) - 20) / (n_randint(state, 2) ? 1.0 : 2.0);
            for (k = 1; k < n; k++)
                y[k] = 0.0;
        }

        eval_plain(z, f, n, x, y);

        arb_init(ax); arb_init(ay); arb_init(az); arb_init(af); arb_init(err); arb_init(tol);
        arf_init(t);
        _dfloat_get_arb(ax, x, n, 0.0);
        _dfloat_get_arb(ay, y, n, 0.0);
        _dfloat_get_arb(az, z, n, 0.0);

        if (arb_is_finite(ax) && arb_is_finite(ay) && eval_arb(af, f, ax, ay, 1500) && arf_cmpabs_2exp_si(arb_midref(af), 800) < 0
            && (arb_is_zero(af) || arf_cmpabs_2exp_si(arb_midref(af), -800) > 0))
        {
            arb_sub(err, az, af, 1500);
            arb_abs(err, err);
            arb_abs(tol, af);
            if (f <= 1)
                arf_set_d(t, f == 0 ? DFLOAT_EXPM1_EPS(n) : DFLOAT_LOG1P_EPS(n));
            else
                arf_set_d(t, ldexp(1.0, 14 - 53 * n));
            arb_mul_arf(tol, tol, t, 1500);
            /* tan near its poles and the others near the zeros of a
               denominator lose relative accuracy as the parts do:
               allow an absolute term of 2^(14 - 53 n) (|x| + 1) */
            if (f == 2)
            {
                arf_set_d(t, ldexp(1.0, 14 - 53 * n) * (fabs(x[0]) + 1.0));
                arb_add_arf(tol, tol, t, 1500);
            }
            if (!arb_is_finite(az) || !arb_le(err, tol))
            {
                flint_printf("FAIL: %s bound (n = %d)\n", fnames[f], n);
                flint_printf("x = %a %a %a %a\ny = %a %a %a %a\n", x[0], n > 1 ? x[1] : 0.0, n > 2 ? x[2] : 0.0, n > 3 ? x[3] : 0.0, y[0], n > 1 ? y[1] : 0.0, n > 2 ? y[2] : 0.0, n > 3 ? y[3] : 0.0);
                flint_printf("x = "); arb_printd(ax, 40); flint_printf("\n");
                flint_printf("y = "); arb_printd(ay, 40); flint_printf("\n");
                flint_printf("z = "); arb_printd(az, 40); flint_printf("\n");
                flint_printf("f = "); arb_printd(af, 40); flint_printf("\n");
                flint_printf("err = "); arb_printd(err, 10); flint_printf(", tol = "); arb_printd(tol, 10); flint_printf("\n");
                flint_abort();
            }
        }

        arb_clear(ax); arb_clear(ay); arb_clear(az); arb_clear(af); arb_clear(err); arb_clear(tol);
        arf_clear(t);
    }

    /* balls: containment at sample points */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        int n, f, k, status;
        double x[DFLOAT_MAX_N + 1], y[DFLOAT_MAX_N + 1], z[DFLOAT_MAX_N + 1];
        arb_t ax, ay, az, af;
        arf_t px, py, t;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        f = n_randint(state, NFUN);
        if (n_randint(state, 5) == 0)
            _dfloat_randtest_special(x, n, state);
        else
            random_arg_el(x, n, f, state);
        random_arg_el(y, n, f == 9 ? 0 : f, state);
        x[n] = n_randint(state, 2) ? 0.0 : (n_randint(state, 2) ? fabs(d_randtest_signed(state, -200, -20)) * (1.0 + fabs(x[0])) : _dfloat_randtest_rad(state));
        y[n] = n_randint(state, 2) ? 0.0 : fabs(d_randtest_signed(state, -200, -20)) * (1.0 + fabs(y[0]));

        status = eval_ball(z, f, n, x, y);

        arb_init(ax); arb_init(ay); arb_init(az); arb_init(af);
        arf_init(px); arf_init(py); arf_init(t);
        _dfloat_get_arb(ax, x, n, x[n]);
        _dfloat_get_arb(ay, y, n, y[n]);
        if (status == GR_SUCCESS)
        {
            _dfloat_get_arb(az, z, n, z[n]);
            for (k = 0; k < 4; k++)
            {
                arf_set(px, arb_midref(ax));
                arf_set(py, arb_midref(ay));
                if (k > 0 && x[n] != D_INF)
                {
                    arf_set_d(t, x[n] * ((k == 1) ? 1.0 : (k == 2) ? -1.0 : (double) n_randint(state, 1000) / 1000.0));
                    arf_add(px, px, t, ARF_PREC_EXACT, ARF_RND_DOWN);
                }
                if (k > 0)
                {
                    arf_set_d(t, y[n] * ((k == 1) ? -1.0 : (k == 2) ? 1.0 : (double) n_randint(state, 1000) / 1000.0));
                    arf_add(py, py, t, ARF_PREC_EXACT, ARF_RND_DOWN);
                }
                if (!arf_is_finite(px) || !arf_is_finite(py))
                    continue;
                arb_set_arf(af, px);
                {
                    arb_t bpy;
                    arb_init(bpy);
                    arb_set_arf(bpy, py);
                    if (!eval_arb(af, f, af, bpy, 2000))
                    {
                        arb_clear(bpy);
                        continue;
                    }
                    arb_clear(bpy);
                }
                if (!arb_contains(az, af))
                {
                    flint_printf("FAIL: %s containment (n = %d, k = %d)\n", fnames[f], n, k);
                    flint_printf("x = "); arb_printd(ax, 30); flint_printf("\n");
                    flint_printf("y = "); arb_printd(ay, 30); flint_printf("\n");
                    flint_printf("z = "); arb_printd(az, 30); flint_printf("\n");
                    flint_printf("f = "); arb_printd(af, 30); flint_printf("\n");
                    flint_abort();
                }
            }
        }
        arb_clear(ax); arb_clear(ay); arb_clear(az); arb_clear(af);
        arf_clear(px); arf_clear(py); arf_clear(t);
    }

    /* the dispatchers of expm1 and log1p (by N, plain and balls) give
       the results of the generic methods */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, ctxb;
        double x[DFLOAT_MAX_N + 1], y[DFLOAT_MAX_N + 1], z[DFLOAT_MAX_N + 1], ry;
        int n = 1 + n_randint(state, DFLOAT_MAX_N), k, st1, st2;

        if (n_randint(state, 8) == 0)
            _dfloat_randtest_special(x, n, state);
        else
            _dfloat_randtest(x, n, state);
        x[n] = n_randint(state, 2) ? 0.0 : ldexp(1.0, -(int) n_randint(state, 200));
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, 0));
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctxb, n, DFLOAT_BALL));

        for (k = 0; k < 4; k++)
        {
            st1 = GR_SUCCESS;
            if (k == 0) { _dfloat_expm1(y, n, x); st2 = gr_expm1(z, x, ctx); }
            else if (k == 1) { _dfloat_log1p(y, n, x); st2 = gr_log1p(z, x, ctx); }
            else if (k == 2) { _dfloat_ball_expm1(y, &ry, n, x, x[n]); y[n] = ry; st2 = gr_expm1(z, x, ctxb); }
            else { st1 = _dfloat_ball_log1p(y, &ry, n, x, x[n]); y[n] = ry; st2 = gr_log1p(z, x, ctxb); }

            /* (the plain methods return GR_DOMAIN where the typed
               functions return nan) */
            if ((k >= 2 && st1 != st2) || (st2 == GR_SUCCESS &&
                memcmp(y, z, sizeof(double) * (n + (k >= 2))) != 0))
            {
                flint_printf("FAIL: expm1 / log1p dispatchers (k = %d, n = %d, status %d %d)\n", k, n, st1, st2);
                flint_printf("x = %a %a %a %a, y = %a %a, z = %a %a\n", x[0], x[1 % n], x[2 % n], x[n], y[0], y[1 % n], z[0], z[1 % n]);
                flint_abort();
            }
        }

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctxb);
    }

    TEST_FUNCTION_END(state);
}
