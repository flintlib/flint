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

/* Rigor test for the ball arithmetic over the whole double range:
   the dNb result must contain the images (computed exactly, or in arb
   at high precision) of sample points of the input balls, including
   the endpoints.  A zero radius must come with an exactly correct
   midpoint, and simple exact inputs must give exact results. */

/* the eight operations under test, on generic (n-double) balls */
typedef struct
{
    double d[DFLOAT_MAX_N + 1];
}
ball_t;

#define RAD(b, n) ((b)->d[n])

static void
ball_randtest(ball_t * b, int n, flint_rand_t state, int special)
{
    switch (n)
    {
        case 1: if (special) d1b_randtest_special((d1b_ptr) b, state); else d1b_randtest((d1b_ptr) b, state); break;
        case 2: if (special) d2b_randtest_special((d2b_ptr) b, state); else d2b_randtest((d2b_ptr) b, state); break;
        case 3: if (special) d3b_randtest_special((d3b_ptr) b, state); else d3b_randtest((d3b_ptr) b, state); break;
        default: if (special) d4b_randtest_special((d4b_ptr) b, state); else d4b_randtest((d4b_ptr) b, state); break;
    }
}

#define OP_ADD 0
#define OP_SUB 1
#define OP_MUL 2
#define OP_SQR 3
#define OP_DIV 4
#define OP_SQRT 5
#define OP_MUL2EXP 6
#define OP_ABS 7
#define OP_MULD 8
#define OP_INV 9
#define OP_RSQRT 10
#define NUM_OPS 11

static const char * op_names[] = { "add", "sub", "mul", "sqr", "div", "sqrt", "mul_2exp", "abs", "mul_d", "inv", "rsqrt" };

#define DISPATCH(n, op, res, x, y, e, dd, status) \
    do { \
        switch (n) { \
        case 1: status = do_op_d1(op, (d1b_ptr) res, (d1b_srcptr) x, (d1b_srcptr) y, e, dd); break; \
        case 2: status = do_op_d2(op, (d2b_ptr) res, (d2b_srcptr) x, (d2b_srcptr) y, e, dd); break; \
        case 3: status = do_op_d3(op, (d3b_ptr) res, (d3b_srcptr) x, (d3b_srcptr) y, e, dd); break; \
        default: status = do_op_d4(op, (d4b_ptr) res, (d4b_srcptr) x, (d4b_srcptr) y, e, dd); break; \
        } \
    } while (0)

#define DEF_DO_OP(X) \
static int do_op_##X(int op, X##b_ptr res, X##b_srcptr x, X##b_srcptr y, slong e, double dd) \
{ \
    switch (op) \
    { \
        case OP_ADD: X##b_add(res, x, y); return GR_SUCCESS; \
        case OP_SUB: X##b_sub(res, x, y); return GR_SUCCESS; \
        case OP_MUL: X##b_mul(res, x, y); return GR_SUCCESS; \
        case OP_SQR: X##b_sqr(res, x); return GR_SUCCESS; \
        case OP_DIV: return X##b_div(res, x, y); \
        case OP_SQRT: return X##b_sqrt(res, x); \
        case OP_MUL2EXP: X##b_mul_2exp_si(res, x, e); return GR_SUCCESS; \
        case OP_ABS: X##b_abs(res, x); return GR_SUCCESS; \
        case OP_MULD: X##b_mul_d(res, x, dd); return GR_SUCCESS; \
        case OP_RSQRT: return X##b_rsqrt(res, x); \
        default: return X##b_inv(res, x); \
    } \
}

DEF_DO_OP(d1)
DEF_DO_OP(d2)
DEF_DO_OP(d3)
DEF_DO_OP(d4)

/* reference: the operation on exact points, in arf when the operation
   is exact and otherwise in arb at high precision */
static int
ref_op_point(arb_t r, int op, const arf_t x, const arf_t y, slong e, double dd, slong prec)
{
    arf_t t;
    arf_init(t);
    switch (op)
    {
        case OP_ADD:
            arf_add(arb_midref(r), x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(r));
            break;
        case OP_SUB:
            arf_sub(arb_midref(r), x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(r));
            break;
        case OP_MUL:
            arf_mul(arb_midref(r), x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(r));
            break;
        case OP_SQR:
            arf_mul(arb_midref(r), x, x, ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(r));
            break;
        case OP_MUL2EXP:
            arf_mul_2exp_si(arb_midref(r), x, e);
            mag_zero(arb_radref(r));
            break;
        case OP_ABS:
            arf_abs(arb_midref(r), x);
            mag_zero(arb_radref(r));
            break;
        case OP_MULD:
            arf_set_d(t, dd);
            arf_mul(arb_midref(r), x, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(r));
            break;
        case OP_DIV:
            arb_set_arf(r, x);
            arb_div_arf(r, r, y, prec);
            break;
        case OP_INV:
            arb_set_arf(r, x);
            arb_inv(r, r, prec);
            break;
        case OP_SQRT:
            arb_sqrt_arf(r, x, prec);
            break;
        case OP_RSQRT:
            arb_set_arf(r, x);
            arb_rsqrt(r, r, prec);
            break;
    }
    arf_clear(t);
    return GR_SUCCESS;
}

/* the status the operation must return, from the exact input balls;
   inflate > 0 widens the radii slightly to detect borderline cases,
   where the implementation may answer GR_UNABLE conservatively */
static int
ref_status(int op, const arb_t x, double rx, const arb_t y, double ry, int inflate)
{
    arf_t lo, hi, r;
    int res = GR_SUCCESS;
    arb_srcptr b = (op == OP_DIV) ? y : x;
    double rad = (op == OP_DIV) ? ry : rx;

    if (op != OP_DIV && op != OP_INV && op != OP_SQRT && op != OP_RSQRT)
        return GR_SUCCESS;

    arf_init(lo); arf_init(hi); arf_init(r);
    if (rad == D_INF)
        arf_pos_inf(r);
    else
    {
        arf_set_d(r, rad);
        if (inflate)
        {
            arf_t t;
            arf_init(t);
            arf_mul_2exp_si(t, r, -40);
            arf_add(r, r, t, ARF_PREC_EXACT, ARF_RND_UP);
            arf_clear(t);
        }
    }
    arf_sub(lo, arb_midref(b), r, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(hi, arb_midref(b), r, ARF_PREC_EXACT, ARF_RND_UP);

    if (op == OP_SQRT)
    {
        if (arf_sgn(hi) < 0) res = GR_DOMAIN;
        else if (arf_sgn(lo) < 0) res = GR_UNABLE;
    }
    else if (op == OP_RSQRT)
    {
        if (arf_sgn(hi) <= 0) res = GR_DOMAIN;
        else if (arf_sgn(lo) <= 0) res = GR_UNABLE;
    }
    else
    {
        if (arf_is_zero(arb_midref(b)) && rad == 0.0) res = GR_DOMAIN;
        else if (arf_sgn(lo) <= 0 && arf_sgn(hi) >= 0) res = GR_UNABLE;
    }
    arf_clear(lo); arf_clear(hi); arf_clear(r);
    return res;
}

/* a point of the ball [mid +/- rad]: the midpoint, an endpoint, or
   something in between; for an infinite radius, any finite point */
static void
sample_point(arf_t p, const arb_t b, double rad, int which, flint_rand_t state)
{
    arf_t t;
    arf_init(t);
    if (rad == D_INF)
    {
        if (which == 0)
            arf_set(p, arb_midref(b));
        else
        {
            arf_set_d(t, d_randtest_signed(state, -1074, 1023));
            arf_add(p, arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            if (!arf_is_finite(p))
                arf_set(p, t);
        }
    }
    else
    {
        arf_set_d(t, rad);
        if (which == 1)
            arf_add(p, arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else if (which == 2)
            arf_sub(p, arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else if (which == 3)
        {
            arf_mul_ui(t, t, n_randint(state, 1000), 30, ARF_RND_DOWN);
            arf_mul_2exp_si(t, t, -10);
            if (n_randint(state, 2))
                arf_neg(t, t);
            arf_add(p, arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        else
            arf_set(p, arb_midref(b));
    }
    arf_clear(t);
}

static void
ball_get_arb(arb_t r, const ball_t * b, int n)
{
    _dfloat_get_arb(r, b->d, n, RAD(b, n));
}

static void
ball_print(const ball_t * b, int n)
{
    int i;
    flint_printf("[");
    for (i = 0; i < n; i++)
        flint_printf("%s%.17g", i ? ", " : "", b->d[i]);
    flint_printf(" +/- %.17g]", RAD(b, n));
}

static int
ball_is_valid(const ball_t * b, int n)
{
    int i;
    for (i = 0; i < n; i++)
        if (!(fabs(b->d[i]) <= DBL_MAX))
            return 0;
    if (!(RAD(b, n) >= 0.0))
        return 0;
    if (b->d[0] == 0.0)
        for (i = 1; i < n; i++)
            if (b->d[i] != 0.0)
                return 0;
    return 1;
}

TEST_FUNCTION_START(ball_arith, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 200000 * flint_test_multiplier(); iter++)
    {
        int n, op, special, status, ref, ref2, aliasing;
        ball_t x, y, r, r2;
        arb_t ax, ay, ar, arr;
        arf_t exact;
        slong e = 0, prec;
        double dd = 0.0;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        op = n_randint(state, NUM_OPS);
        special = (n_randint(state, 4) == 0);
        aliasing = n_randint(state, 4);
        prec = 2000 + n_randint(state, 100);

        ball_randtest(&x, n, state, special);
        ball_randtest(&y, n, state, special);
        if (n_randint(state, 8) == 0)
            y = x;   /* correlated inputs: cancellation */
        if (op == OP_MUL2EXP)
            e = (slong) n_randint(state, 4200) - 2100;
        if (op == OP_MULD)
            dd = special ? d_randtest_signed(state, -1074, 1023) : d_randtest_signed(state, -20, 20);

        arb_init(ax); arb_init(ay); arb_init(ar); arb_init(arr);
        arf_init(exact);
        ball_get_arb(ax, &x, n);
        ball_get_arb(ay, &y, n);

        ref = ref_status(op, ax, RAD(&x, n), ay, RAD(&y, n), 0);
        ref2 = ref_status(op, ax, RAD(&x, n), ay, RAD(&y, n), 1);

        /* run the operation, possibly with aliasing */
        if (aliasing == 1)
        {
            r = x;
            DISPATCH(n, op, &r, &r, &y, e, dd, status);
        }
        else if (aliasing == 2)
        {
            r = y;
            DISPATCH(n, op, &r, &x, &r, e, dd, status);
        }
        else
        {
            DISPATCH(n, op, &r, &x, &y, e, dd, status);
        }
        /* the non-aliased result must agree bitwise */
        DISPATCH(n, op, &r2, &x, &y, e, dd, status);
        if (status == GR_SUCCESS && memcmp(&r, &r2, sizeof(double) * (n + 1)) != 0)
        {
            flint_printf("FAIL: aliasing\n");
            flint_printf("n = %d, op = %s\n", n, op_names[op]);
            ball_print(&x, n); flint_printf("\n");
            ball_print(&y, n); flint_printf("\n");
            ball_print(&r, n); flint_printf("\n");
            ball_print(&r2, n); flint_printf("\n");
            flint_abort();
        }

        /* GR_UNABLE is acceptable in borderline cases */
        if (status != ref && !(status == GR_UNABLE && ref != ref2))
        {
            flint_printf("FAIL: status\n");
            flint_printf("n = %d, op = %s, status = %d, ref = %d, e = %wd\n", n, op_names[op], status, ref, e);
            ball_print(&x, n); flint_printf("\n");
            ball_print(&y, n); flint_printf("\n");
            flint_abort();
        }

        if (status == GR_SUCCESS)
        {
            int kx, ky;
            arf_t px, py;

            if (!ball_is_valid(&r, n))
            {
                flint_printf("FAIL: invalid result\n");
                flint_printf("n = %d, op = %s\n", n, op_names[op]);
                ball_print(&x, n); flint_printf("\n");
                ball_print(&y, n); flint_printf("\n");
                ball_print(&r, n); flint_printf("\n");
                flint_abort();
            }

            ball_get_arb(arr, &r, n);
            arf_init(px); arf_init(py);

            /* containment of the images of sample points */
            for (kx = 0; kx < 4; kx++)
            for (ky = 0; ky < 4; ky++)
            {
                sample_point(px, ax, RAD(&x, n), kx, state);
                sample_point(py, ay, RAD(&y, n), ky, state);

                /* the point must lie in the domain */
                if (op == OP_DIV && arf_is_zero(py)) continue;
                if (op == OP_INV && arf_is_zero(px)) continue;
                if (op == OP_SQRT && arf_sgn(px) < 0) continue;
                if (op == OP_RSQRT && arf_sgn(px) <= 0) continue;

                ref_op_point(ar, op, px, py, e, dd, prec);

                if (!arb_contains(arr, ar))
                {
                    flint_printf("FAIL: containment\n");
                    flint_printf("n = %d, op = %s, e = %wd, dd = %.17g, kx = %d, ky = %d\n", n, op_names[op], e, dd, kx, ky);
                    flint_printf("x = "); ball_print(&x, n); flint_printf("\n");
                    flint_printf("y = "); ball_print(&y, n); flint_printf("\n");
                    flint_printf("r = "); ball_print(&r, n); flint_printf("\n");
                    flint_printf("px = "); arf_printd(px, 50); flint_printf("\n");
                    flint_printf("py = "); arf_printd(py, 50); flint_printf("\n");
                    flint_printf("ar = "); arb_printd(ar, 50); flint_printf("\n");
                    flint_printf("arr = "); arb_printd(arr, 50); flint_printf("\n");
                    flint_abort();
                }

                /* an exact result must be the exact image of the midpoints */
                if (kx == 0 && ky == 0 && RAD(&r, n) == 0.0 && arb_is_exact(ar)
                    && !arf_equal(arb_midref(arr), arb_midref(ar)))
                {
                    flint_printf("FAIL: claimed exact but wrong\n");
                    flint_printf("n = %d, op = %s\n", n, op_names[op]);
                    flint_printf("x = "); ball_print(&x, n); flint_printf("\n");
                    flint_printf("y = "); ball_print(&y, n); flint_printf("\n");
                    flint_printf("r = "); ball_print(&r, n); flint_printf("\n");
                    flint_printf("exact = "); arf_printd(arb_midref(ar), 50); flint_printf("\n");
                    flint_abort();
                }
            }

            /* a finite exact result with few enough bits from short
               exact inputs in a comfortable range must come out exact */
            if (arb_is_exact(ax) && arb_is_exact(ay)
                && (op == OP_ADD || op == OP_SUB || op == OP_MUL || op == OP_SQR || op == OP_MULD))
            {
                ref_op_point(ar, op, arb_midref(ax), arb_midref(ay), e, dd, prec);
                arf_set(exact, arb_midref(ar));
                if (arf_is_finite(exact) && !arf_is_zero(exact)
                    && arf_bits(exact) <= 40 && arf_cmpabs_2exp_si(exact, 400) < 0
                    && arf_cmpabs_2exp_si(exact, -400) > 0
                    && arf_bits(arb_midref(ax)) <= 40 && arf_bits(arb_midref(ay)) <= 40
                    && arf_cmpabs_2exp_si(arb_midref(ax), 300) < 0 && arf_cmpabs_2exp_si(arb_midref(ay), 300) < 0
                    && (arf_is_zero(arb_midref(ax)) || arf_cmpabs_2exp_si(arb_midref(ax), -300) > 0)
                    && (arf_is_zero(arb_midref(ay)) || arf_cmpabs_2exp_si(arb_midref(ay), -300) > 0)
                    && (op != OP_MULD || (fabs(dd) < 1e30 && (fabs(dd) > 1e-30 || dd == 0.0)))
                    && RAD(&r, n) != 0.0)
                {
                    flint_printf("FAIL: expected an exact result\n");
                    flint_printf("n = %d, op = %s\n", n, op_names[op]);
                    flint_printf("x = "); ball_print(&x, n); flint_printf("\n");
                    flint_printf("y = "); ball_print(&y, n); flint_printf("\n");
                    flint_printf("r = "); ball_print(&r, n); flint_printf("\n");
                    flint_abort();
                }
            }

            arf_clear(px); arf_clear(py);
        }

        arb_clear(ax); arb_clear(ay); arb_clear(ar); arb_clear(arr);
        arf_clear(exact);
    }

    TEST_FUNCTION_END(state);
}
