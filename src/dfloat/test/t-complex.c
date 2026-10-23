/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"

/* The complex types dNc and dNcb: the generic ring axioms; the ball
   functions against acb (containment of the exact value at the
   midpoint, and the statuses on and off the domains); the plain
   functions against acb (normwise relative accuracy); the vector
   operations and dot products against scalar loops; conversions,
   strings and the special cases of the conventions. */

#define C_CHECK(cond, msg) \
    do { \
        if (!(cond)) \
        { \
            flint_printf("FAIL: %s (line %d), n = %d, flags = %d\n", msg, __LINE__, n, flags); \
            flint_abort(); \
        } \
    } while (0)

#define NFUNC 16

/* f(x) (and f(x, y) for pow, div) in the generic ring and in acb */
static int
_t_c_eval(int f, gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    switch (f)
    {
        case 0: return gr_exp(res, x, ctx);
        case 1: return gr_log(res, x, ctx);
        case 2: return gr_sqrt(res, x, ctx);
        case 3: return gr_rsqrt(res, x, ctx);
        case 4: return gr_sin(res, x, ctx);
        case 5: return gr_cos(res, x, ctx);
        case 6: return gr_tan(res, x, ctx);
        case 7: return gr_sinh(res, x, ctx);
        case 8: return gr_cosh(res, x, ctx);
        case 9: return gr_tanh(res, x, ctx);
        case 10: return gr_inv(res, x, ctx);
        case 11: return gr_div(res, x, y, ctx);
        case 12: return gr_pow(res, x, y, ctx);
        case 13: return gr_mul(res, x, y, ctx);
        case 14: return gr_sqr(res, x, ctx);
        default: return gr_abs(res, x, ctx);
    }
}

static void
_t_c_ref(int f, acb_t res, const acb_t x, const acb_t y, slong prec)
{
    switch (f)
    {
        case 0: acb_exp(res, x, prec); break;
        case 1: acb_log(res, x, prec); break;
        case 2: acb_sqrt(res, x, prec); break;
        case 3: acb_rsqrt(res, x, prec); break;
        case 4: acb_sin(res, x, prec); break;
        case 5: acb_cos(res, x, prec); break;
        case 6: acb_tan(res, x, prec); break;
        case 7: acb_sinh(res, x, prec); break;
        case 8: acb_cosh(res, x, prec); break;
        case 9: acb_tanh(res, x, prec); break;
        case 10: acb_inv(res, x, prec); break;
        case 11: acb_div(res, x, y, prec); break;
        case 12: acb_pow(res, x, y, prec); break;
        case 13: acb_mul(res, x, y, prec); break;
        case 14: acb_sqr(res, x, prec); break;
        default: acb_abs(acb_realref(res), x, prec); arb_zero(acb_imagref(res)); break;
    }
}

/* a random complex number with parts of moderate size (sometimes zero,
   sometimes a small integer) */
static void
_t_c_rand(acb_t x, flint_rand_t state, int mag)
{
    int k;
    for (k = 0; k < 2; k++)
    {
        arb_ptr p = k ? acb_imagref(x) : acb_realref(x);
        switch (n_randint(state, 6))
        {
            case 0: arb_zero(p); break;
            case 1: arb_set_si(p, (slong) n_randint(state, 7) - 3); break;
            default:
                arf_set_d(arb_midref(p), ldexp((double) n_randlimb(state) * 0x1p-64 - 0.5,
                    (int) n_randint(state, 2 * mag + 1) - mag + 1));
                mag_zero(arb_radref(p));
        }
    }
}

/* the dfloat ball of the exact value x plus a random radius */
static void
_t_c_ball(gr_ptr res, const acb_t x, gr_ctx_t ctx, gr_ctx_t CC, flint_rand_t state)
{
    acb_t t;
    acb_init(t);
    acb_set(t, x);
    if (n_randint(state, 2))
    {
        mag_set_d(arb_radref(acb_realref(t)), ldexp(1.0, -(int) n_randint(state, 60)));
        mag_set_d(arb_radref(acb_imagref(t)), ldexp(1.0, -(int) n_randint(state, 60)));
    }
    GR_MUST_SUCCEED(gr_set_other(res, t, CC, ctx));
    acb_clear(t);
}

/* the rings' conversions of an element to acb */
static void
_t_c_get(acb_t res, gr_srcptr x, gr_ctx_t ctx, gr_ctx_t CC)
{
    GR_MUST_SUCCEED(gr_set_other(res, x, ctx, CC));
}

static void
_t_c_functions(gr_ctx_t ctx, int n, int flags, flint_rand_t state, slong iter)
{
    gr_ctx_t CC;
    gr_ptr x, y, z;
    acb_t a, b, r, s;
    slong i;
    int f, status, ball = (flags & DFLOAT_BALL) != 0;

    gr_ctx_init_complex_acb(CC, 400);
    x = gr_heap_init(ctx); y = gr_heap_init(ctx); z = gr_heap_init(ctx);
    acb_init(a); acb_init(b); acb_init(r); acb_init(s);

    for (i = 0; i < iter; i++)
    {
        f = (int) n_randint(state, NFUNC);
        _t_c_rand(a, state, (f == 0 || f == 4 || f == 5 || f == 7 || f == 8) ? 5 : 20);
        _t_c_rand(b, state, (f == 12) ? 3 : 20);
        if (ball)
        {
            _t_c_ball(x, a, ctx, CC, state);
            _t_c_ball(y, b, ctx, CC, state);
            /* the exact points: the midpoints */
            _t_c_get(a, x, ctx, CC);
            _t_c_get(b, y, ctx, CC);
            arb_get_mid_arb(acb_realref(a), acb_realref(a));
            arb_get_mid_arb(acb_imagref(a), acb_imagref(a));
            arb_get_mid_arb(acb_realref(b), acb_realref(b));
            arb_get_mid_arb(acb_imagref(b), acb_imagref(b));
        }
        else
        {
            GR_MUST_SUCCEED(gr_set_other(x, a, CC, ctx));
            GR_MUST_SUCCEED(gr_set_other(y, b, CC, ctx));
            _t_c_get(a, x, ctx, CC);
            _t_c_get(b, y, ctx, CC);
        }

        status = _t_c_eval(f, z, x, y, ctx);
        _t_c_ref(f, r, a, b, 400);

        if (ball)
        {
            if (status == GR_SUCCESS)
            {
                _t_c_get(s, z, ctx, CC);
                /* (the reference with more precision where it lost some) */
                if (acb_is_finite(r) && !acb_contains(s, r))
                    _t_c_ref(f, r, a, b, 4000);
                if (acb_is_finite(r) && !acb_contains(s, r))
                {
                    flint_printf("FAIL: containment, f = %d\n", f);
                    flint_printf("x = "); gr_println(x, ctx);
                    flint_printf("y = "); gr_println(y, ctx);
                    flint_printf("z = "); gr_println(z, ctx);
                    flint_printf("r = "); acb_printn(r, 30, 0); flint_printf("\n");
                    flint_abort();
                }
            }
            /* the domains: x = 0 for log, rsqrt, inv; y = 0 for div */
            if ((f == 1 || f == 3 || f == 10) && acb_is_zero(a) && gr_is_zero(x, ctx) == T_TRUE)
                C_CHECK(status == GR_DOMAIN, "domain");
            else if (f == 11 && gr_is_zero(y, ctx) == T_TRUE)
                C_CHECK(status == GR_DOMAIN, "domain");
            else if (status == GR_DOMAIN)
                C_CHECK(f == 12 && acb_is_zero(a), "unexpected GR_DOMAIN");
            if (f == 0 || f == 2 || f == 4 || f == 5 || f == 7 || f == 8 || f == 13 || f == 14 || f == 15)
                C_CHECK(status == GR_SUCCESS, "entire");
        }
        else if (status == GR_SUCCESS && acb_is_finite(r))
        {
            /* normwise: |z - f| <= 2^(24 - 53 n) (1 + |y log x|) max(|f|, tiny) */
            arb_t e, m;
            arb_init(e); arb_init(m);
            if (gr_set_other(s, z, ctx, CC) != GR_SUCCESS)
            {
                flint_printf("FAIL: nonfinite result, f = %d\n", f);
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("z = "); gr_println(z, ctx);
                flint_printf("r = "); acb_printn(r, 30, 0); flint_printf("\n");
                flint_abort();
            }
            acb_sub(s, s, r, 400);
            acb_abs(e, s, 400);
            acb_abs(m, r, 400);
            if (f == 12)
            {
                acb_t l;
                acb_init(l);
                acb_log(l, a, 400);
                acb_mul(l, l, b, 400);
                acb_abs(acb_realref(s), l, 400);
                arb_add_ui(acb_realref(s), acb_realref(s), 1, 400);
                arb_mul(m, m, acb_realref(s), 400);
                acb_clear(l);
            }
            arb_mul_2exp_si(m, m, 24 - 53 * n);
            if (f != 6 && f != 9 && arf_cmpabs_2exp_si(arb_midref(m), -900) > 0 && arb_gt(e, m))
            {
                flint_printf("FAIL: accuracy, f = %d\n", f);
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("z = "); gr_println(z, ctx);
                flint_printf("r = "); acb_printn(r, 30, 0); flint_printf("\n");
                flint_abort();
            }
            arb_clear(e); arb_clear(m);
        }
    }

    gr_heap_clear(x, ctx); gr_heap_clear(y, ctx); gr_heap_clear(z, ctx);
    acb_clear(a); acb_clear(b); acb_clear(r); acb_clear(s);
    gr_ctx_clear(CC);
}

/* the vector operations and dot products against scalar loops (the
   scalar results must be contained in (balls) or close to (plain) the
   vector ones) */
static void
_t_c_vectors(gr_ctx_t ctx, int n, int flags, flint_rand_t state, slong iter)
{
    gr_ctx_t CC;
    gr_ptr x, y, z, w, c, t, u;
    acb_t a, b;
    slong i, j, len;
    int op, ball = (flags & DFLOAT_BALL) != 0, sub;
    slong sz = ctx->sizeof_elem;

    gr_ctx_init_complex_acb(CC, 400);
    acb_init(a); acb_init(b);
    x = gr_heap_init_vec(80, ctx); y = gr_heap_init_vec(80, ctx);
    z = gr_heap_init_vec(80, ctx); w = gr_heap_init_vec(80, ctx);
    c = gr_heap_init(ctx); t = gr_heap_init(ctx); u = gr_heap_init(ctx);

    for (i = 0; i < iter; i++)
    {
        len = n_randint(state, 80);
        op = (int) n_randint(state, 7);
        sub = (int) n_randint(state, 2);
        for (j = 0; j < len; j++)
        {
            _t_c_rand(a, state, (op == 6) ? 5 : 20);
            if (ball) _t_c_ball(GR_ENTRY(x, j, sz), a, ctx, CC, state);
            else GR_MUST_SUCCEED(gr_set_other(GR_ENTRY(x, j, sz), a, CC, ctx));
            _t_c_rand(a, state, 20);
            if (ball) _t_c_ball(GR_ENTRY(y, j, sz), a, ctx, CC, state);
            else GR_MUST_SUCCEED(gr_set_other(GR_ENTRY(y, j, sz), a, CC, ctx));
        }
        _t_c_rand(a, state, 20);
        GR_MUST_SUCCEED(gr_set_other(c, a, CC, ctx));
        GR_MUST_SUCCEED(_gr_vec_set(z, y, len, ctx));
        GR_MUST_SUCCEED(_gr_vec_set(w, y, len, ctx));

        /* z: the vector operation, w: the scalar loop */
        switch (op)
        {
            case 0: GR_MUST_SUCCEED(_gr_vec_add(z, x, y, len, ctx)); break;
            case 1: GR_MUST_SUCCEED(_gr_vec_mul(z, x, y, len, ctx)); break;
            case 2: GR_MUST_SUCCEED(_gr_vec_mul_scalar(z, x, len, c, ctx)); break;
            case 3: GR_MUST_SUCCEED(_gr_vec_addmul_scalar(z, x, len, c, ctx)); break;
            case 4: GR_MUST_SUCCEED(_gr_vec_submul_scalar(z, x, len, c, ctx)); break;
            case 5: GR_MUST_SUCCEED(sub ? _gr_vec_dot_rev(z, c, 1, x, y, len, ctx)
                                        : _gr_vec_dot(z, c, 0, x, y, len, ctx)); break;
            default: GR_MUST_SUCCEED(_gr_vec_exp(z, x, len, ctx)); break;
        }
        if (op == 5)
        {
            GR_MUST_SUCCEED(gr_set(t, c, ctx));
            for (j = 0; j < len; j++)
            {
                GR_MUST_SUCCEED(gr_mul(u, GR_ENTRY(x, j, sz), GR_ENTRY(y, sub ? len - 1 - j : j, sz), ctx));
                GR_MUST_SUCCEED(sub ? gr_sub(t, t, u, ctx) : gr_add(t, t, u, ctx));
            }
            GR_MUST_SUCCEED(gr_set(w, t, ctx));
        }
        else
            for (j = 0; j < len; j++)
            {
                gr_ptr xj = GR_ENTRY(x, j, sz), yj = GR_ENTRY(y, j, sz), wj = GR_ENTRY(w, j, sz);
                switch (op)
                {
                    case 0: GR_MUST_SUCCEED(gr_add(wj, xj, yj, ctx)); break;
                    case 1: GR_MUST_SUCCEED(gr_mul(wj, xj, yj, ctx)); break;
                    case 2: GR_MUST_SUCCEED(gr_mul(wj, xj, c, ctx)); break;
                    case 3: GR_MUST_SUCCEED(gr_mul(u, xj, c, ctx)); GR_MUST_SUCCEED(gr_add(wj, yj, u, ctx)); break;
                    case 4: GR_MUST_SUCCEED(gr_mul(u, xj, c, ctx)); GR_MUST_SUCCEED(gr_sub(wj, yj, u, ctx)); break;
                    default: GR_MUST_SUCCEED(gr_exp(wj, xj, ctx)); break;
                }
            }

        for (j = 0; j < ((op == 5) ? 1 : len); j++)
        {
            _t_c_get(a, GR_ENTRY(z, j, sz), ctx, CC);
            _t_c_get(b, GR_ENTRY(w, j, sz), ctx, CC);
            if (ball)
            {
                /* both enclose the exact value at the midpoints */
                C_CHECK(acb_overlaps(a, b), "vector overlap");
            }
            else
            {
                arb_t e, m;
                arb_init(e); arb_init(m);
                acb_sub(a, a, b, 400);
                acb_abs(e, a, 400);
                acb_abs(m, b, 400);
                arb_add_ui(m, m, 1, 400);
                arb_mul_2exp_si(m, m, 30 + (op == 5 ? 10 : 0) - 53 * n);
                if (acb_is_finite(b) && arb_gt(e, m))
                {
                    flint_printf("FAIL: vector op %d, len %wd, index %wd\n", op, len, j);
                    gr_println(GR_ENTRY(z, j, sz), ctx); gr_println(GR_ENTRY(w, j, sz), ctx);
                    flint_abort();
                }
                arb_clear(e); arb_clear(m);
            }
        }
    }

    gr_heap_clear_vec(x, 80, ctx); gr_heap_clear_vec(y, 80, ctx);
    gr_heap_clear_vec(z, 80, ctx); gr_heap_clear_vec(w, 80, ctx);
    gr_heap_clear(c, ctx); gr_heap_clear(t, ctx); gr_heap_clear(u, ctx);
    acb_clear(a); acb_clear(b);
    gr_ctx_clear(CC);
}

/* elements from strings: re, im and the radii (ball) */
static void
_t_c_set(gr_ptr res, gr_ctx_t ctx, double re, double rr, double im, double ir)
{
    double * d = res;
    int i, n = DFLOAT_CTX_N(ctx), ball = DFLOAT_CTX_BALL(ctx);
    for (i = 0; i < 2 * (n + ball); i++)
        d[i] = 0.0;
    d[0] = re;
    d[n + ball] = im;
    if (ball)
    {
        d[n] = rr;
        d[2 * n + 1] = ir;
    }
}

/* the conventions: domains, the whole line in a part, branch cuts,
   exact results, conversions and strings */
static void
_t_c_special(int n, int flags)
{
    gr_ctx_t ctx, R, RR, CC, D;
    gr_ptr x, y, z;
    acb_t a;
    arb_t t;
    fmpz_t f;
    char * s;
    int ball = (flags & DFLOAT_BALL) != 0;
    double * d;

    GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, flags));
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(R, n, flags & ~DFLOAT_COMPLEX));
    gr_ctx_init_real_arb(RR, 300);
    gr_ctx_init_complex_acb(CC, 300);
    x = gr_heap_init(ctx); y = gr_heap_init(ctx); z = gr_heap_init(ctx);
    d = z;
    acb_init(a); arb_init(t); fmpz_init(f);

    /* exact arithmetic: (1 + 2i)(3 - i) = 5 + 5i, (5 + 5i) / (1 + 2i) = 3 - i */
    _t_c_set(x, ctx, 1.0, 0.0, 2.0, 0.0);
    _t_c_set(y, ctx, 3.0, 0.0, -1.0, 0.0);
    GR_MUST_SUCCEED(gr_mul(z, x, y, ctx));
    C_CHECK(d[0] == 5.0 && d[n + ball] == 5.0 && (!ball || (d[n] == 0.0 && d[2 * n + 1] == 0.0)), "exact mul");
    GR_MUST_SUCCEED(gr_div(z, z, x, ctx));
    C_CHECK(gr_equal(z, y, ctx) == T_TRUE || (!ball && fabs(d[0] - 3.0) < 1e-15), "div");
    /* sqrt(-4) = 2i, sqrt(2i) = 1 + i, i^2 = -1 */
    _t_c_set(x, ctx, -4.0, 0.0, 0.0, 0.0);
    GR_MUST_SUCCEED(gr_sqrt(z, x, ctx));
    C_CHECK(d[0] == 0.0 && d[n + ball] == 2.0 && (!ball || d[2 * n + 1] == 0.0), "sqrt(-4)");
    _t_c_set(x, ctx, 0.0, 0.0, 2.0, 0.0);
    GR_MUST_SUCCEED(gr_sqrt(z, x, ctx));
    C_CHECK(d[0] == 1.0 && d[n + ball] == 1.0, "sqrt(2i)");
    GR_MUST_SUCCEED(gr_i(x, ctx));
    GR_MUST_SUCCEED(gr_sqr(z, x, ctx));
    C_CHECK(gr_is_neg_one(z, ctx) == T_TRUE, "i^2");
    /* log(-1) = pi i */
    _t_c_set(x, ctx, -1.0, 0.0, 0.0, 0.0);
    GR_MUST_SUCCEED(gr_log(z, x, ctx));
    C_CHECK(d[0] == 0.0 && fabs(d[n + ball] - 3.141592653589793) < 1e-15, "log(-1)");
    /* exp(0) = 1 exactly; exp(pi i/2) ~ i */
    GR_MUST_SUCCEED(gr_zero(x, ctx));
    GR_MUST_SUCCEED(gr_exp(z, x, ctx));
    C_CHECK(gr_is_one(z, ctx) == T_TRUE, "exp(0)");

    /* domains and statuses */
    GR_MUST_SUCCEED(gr_zero(y, ctx));
    GR_MUST_SUCCEED(gr_one(x, ctx));
    C_CHECK(gr_div(z, x, y, ctx) == GR_DOMAIN, "1/0");
    C_CHECK(gr_inv(z, y, ctx) == GR_DOMAIN, "inv(0)");
    C_CHECK(gr_rsqrt(z, y, ctx) == GR_DOMAIN, "rsqrt(0)");
    _t_c_set(y, ctx, 0.0, 0.0, 2.0, 0.0);
    C_CHECK(gr_pow(z, x, y, ctx) == GR_SUCCESS, "1^(2i)");
    GR_MUST_SUCCEED(gr_zero(x, ctx));
    _t_c_set(y, ctx, 2.0, 0.0, 1.0, 0.0);
    C_CHECK(gr_pow(z, x, y, ctx) == GR_SUCCESS && gr_is_zero(z, ctx) == T_TRUE, "0^(2+i)");
    _t_c_set(y, ctx, -1.0, 0.0, 0.0, 0.0);
    C_CHECK(gr_pow(z, x, y, ctx) == GR_DOMAIN, "0^-1");
    _t_c_set(y, ctx, 0.0, 0.0, 1.0, 0.0);
    C_CHECK(gr_pow(z, x, y, ctx) == GR_DOMAIN, "0^i");
    _t_c_set(y, ctx, 0.0, 0.0, 0.0, 0.0);
    C_CHECK(gr_pow(z, x, y, ctx) == GR_SUCCESS && gr_is_one(z, ctx) == T_TRUE, "0^0");

    if (ball)
    {
        C_CHECK(gr_log(z, x, ctx) == GR_DOMAIN, "log(0)");
        /* balls containing zero */
        _t_c_set(y, ctx, 0.0, 0.5, 0.0, 0.5);
        GR_MUST_SUCCEED(gr_one(x, ctx));
        C_CHECK(gr_div(z, x, y, ctx) == GR_UNABLE, "1/[0 +/- 1/2]");
        C_CHECK(gr_log(z, y, ctx) == GR_UNABLE, "log");
        C_CHECK(gr_rsqrt(z, y, ctx) == GR_UNABLE, "rsqrt");
        C_CHECK(gr_sqrt(z, y, ctx) == GR_SUCCESS, "sqrt (entire)");
        C_CHECK(gr_exp(z, y, ctx) == GR_SUCCESS, "exp (entire)");
        _t_c_set(x, ctx, 0.0, 0.0, 0.0, 0.0);
        _t_c_set(y, ctx, 0.5, 0.25, 0.0, 0.0);
        C_CHECK(gr_pow(z, y, x, ctx) == GR_SUCCESS && gr_is_one(z, ctx) == T_TRUE, "y^0");
        _t_c_set(x, ctx, 0.0, 0.5, 0.0, 0.5);
        C_CHECK(gr_pow(z, x, y, ctx) == GR_SUCCESS, "[0 +/- 1/2]^[1/2 +/- 1/4] (Re y > 0)");
        /* a ball crossing the branch cut: an enclosure of both sides */
        _t_c_set(x, ctx, -1.0, 0.0, 0.0, 0.25);
        GR_MUST_SUCCEED(gr_log(z, x, ctx));
        C_CHECK(d[n + ball] == 0.0 && d[2 * n + 1] >= 3.14159, "log across the cut");
        GR_MUST_SUCCEED(gr_sqrt(z, x, ctx));
        C_CHECK(d[n + ball] == 0.0 && d[2 * n + 1] >= 1.0, "sqrt across the cut");
        /* tan near a real pole off the real line (the pole is not in the ball) */
        _t_c_set(x, ctx, 1.5707963267948966, 0.01, 1e-3, 0.0);
        C_CHECK(gr_tan(z, x, ctx) == GR_SUCCESS, "tan near a pole");
        _t_c_set(x, ctx, 1.5707963267948966, 0.01, 0.0, 1e-3);
        C_CHECK(gr_tan(z, x, ctx) == GR_UNABLE, "tan at a pole");

        /* the whole line in a part */
        _t_c_set(x, ctx, 0.0, D_INF, 0.0, 0.0);          /* [0 +/- inf] + 0i */
        GR_MUST_SUCCEED(gr_set_si(y, 3, ctx));
        C_CHECK(gr_div(z, x, y, ctx) == GR_SUCCESS && d[n] == D_INF && d[2 * n + 1] == 0.0, "W/3");
        C_CHECK(gr_exp(z, x, ctx) == GR_SUCCESS && d[n] == D_INF && gr_is_zero(z, ctx) != T_TRUE, "exp(W)");
        C_CHECK(gr_sin(z, x, ctx) == GR_SUCCESS && d[n] <= 1.0 + 1e-12, "sin(W) in [-1, 1]");
        C_CHECK(gr_div(z, y, x, ctx) == GR_UNABLE, "3/W");
        C_CHECK(gr_log(z, x, ctx) == GR_UNABLE, "log(W)");
        GR_MUST_SUCCEED(gr_zero(y, ctx));
        C_CHECK(gr_div(z, x, y, ctx) == GR_DOMAIN, "W/0");
        C_CHECK(gr_mul(z, x, y, ctx) == GR_SUCCESS && gr_is_zero(z, ctx) == T_TRUE, "W * 0");
        _t_c_set(x, ctx, D_NAN, 0.0, 2.0, 0.0);          /* W + 2i */
        C_CHECK(gr_inv(z, x, ctx) == GR_SUCCESS && d[n] <= 0.5 + 1e-12 && d[2 * n + 1] <= 0.5 + 1e-12, "1/(W + 2i)");
        C_CHECK(gr_log(z, x, ctx) == GR_SUCCESS && d[n] == D_INF, "log(W + 2i)");
        C_CHECK(gr_tan(z, x, ctx) == GR_SUCCESS && d[n] <= 2.0, "tan(W + 2i)");
        C_CHECK(gr_exp(z, x, ctx) == GR_SUCCESS, "exp(W + 2i)");
        GR_MUST_SUCCEED(gr_set(z, x, ctx));
        C_CHECK(d[0] == 0.0 && d[n] == D_INF, "set normalizes");
        C_CHECK(gr_is_zero(x, ctx) == T_FALSE, "W + 2i != 0");
        C_CHECK(gr_equal(x, x, ctx) == T_UNKNOWN, "equal");
        _t_c_set(x, ctx, 0.0, 0.0, D_INF, 0.0);          /* i W: 0 + [0 +/- inf] i */
        C_CHECK(gr_is_zero(x, ctx) == T_UNKNOWN, "is_zero(i W)");
        C_CHECK(gr_exp(z, x, ctx) == GR_SUCCESS && d[n] <= 1.0 + 1e-12 && d[2 * n + 1] <= 1.0 + 1e-12, "exp(i W)");

        /* conversions: real parts from/to the real rings */
        _t_c_set(x, ctx, 2.0, 0.0, 0.0, 0.0);
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_SUCCESS && fmpz_equal_si(f, 2), "get_fmpz");
        _t_c_set(x, ctx, 2.0, 0.0, 1.0, 0.0);
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_DOMAIN, "get_fmpz (complex)");
        C_CHECK(gr_set_other(t, x, ctx, RR) == GR_DOMAIN, "to arb (complex)");
        C_CHECK(gr_set_other(t, x, ctx, R) == GR_DOMAIN, "to real (complex)");
        _t_c_set(x, ctx, 2.0, 0.0, 0.0, 1.0);
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_UNABLE, "get_fmpz (maybe real)");
        C_CHECK(gr_set_other(t, x, ctx, RR) == GR_UNABLE, "to arb (maybe real)");
        _t_c_set(x, ctx, 2.5, 0.0, 0.0, 0.0);
        C_CHECK(gr_set_other(t, x, ctx, RR) == GR_SUCCESS && arb_equal_si(t, 0) == 0, "to arb");
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_DOMAIN, "get_fmpz (not an integer)");
    }
    else
    {
        _t_c_set(x, ctx, 2.0, 0.0, 1.0, 0.0);
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_DOMAIN, "get_fmpz (complex)");
        C_CHECK(gr_set_other(t, x, ctx, RR) == GR_DOMAIN, "to arb (complex)");
        _t_c_set(x, ctx, -3.0, 0.0, 0.0, 0.0);
        C_CHECK(gr_get_fmpz(f, x, ctx) == GR_SUCCESS && fmpz_equal_si(f, -3), "get_fmpz");
    }

    /* conversions between the complex formats, and to acb */
    _t_c_set(x, ctx, 0.1, ball ? 1e-20 : 0.0, -0.3, 0.0);
    {
        int m;
        for (m = 1; m <= DFLOAT_MAX_N; m++)
        {
            gr_ptr v;
            acb_t b;
            GR_MUST_SUCCEED(gr_ctx_init_dfloat(D, m, flags));
            v = gr_heap_init(D);
            acb_init(b);
            GR_MUST_SUCCEED(gr_set_other(v, x, ctx, D));
            GR_MUST_SUCCEED(gr_set_other(a, x, ctx, CC));
            GR_MUST_SUCCEED(gr_set_other(b, v, D, CC));
            if (ball)
                C_CHECK(acb_overlaps(a, b), "conversion between formats");
            else
                C_CHECK(fabs(arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR) - 0.1) < 1e-15, "conversion between formats");
            acb_clear(b);
            gr_heap_clear(v, D);
            gr_ctx_clear(D);
        }
    }
    /* from and to the real ring */
    GR_MUST_SUCCEED(gr_set_si(y, -7, R));
    GR_MUST_SUCCEED(gr_set_other(z, y, R, ctx));
    C_CHECK(gr_is_neg_one(z, ctx) == T_FALSE && d[0] == -7.0 && d[n + ball] == 0.0, "from real");
    GR_MUST_SUCCEED(gr_set_other(y, z, ctx, R));

    /* strings */
    _t_c_set(x, ctx, 1.5, 0.0, -2.5, 0.0);
    GR_MUST_SUCCEED(gr_get_str(&s, x, ctx));
    C_CHECK(strstr(s, "1.5") != NULL && strstr(s, " - 2.5") != NULL && strstr(s, "*I") != NULL, "string");
    flint_free(s);
    GR_MUST_SUCCEED(gr_set_str(y, "1.5 - 2.5*I", ctx));
    C_CHECK(gr_equal(x, y, ctx) == T_TRUE, "string round trip");
    _t_c_set(x, ctx, 0.0, 0.0, 3.0, 0.0);
    GR_MUST_SUCCEED(gr_get_str(&s, x, ctx));
    C_CHECK(strcmp(s, "3.0000000000000000*I") == 0 || strstr(s, "*I") != NULL, "string (imaginary)");
    flint_free(s);
    GR_MUST_SUCCEED(gr_one(x, ctx));
    GR_MUST_SUCCEED(gr_get_str(&s, x, ctx));
    C_CHECK(ball ? (strcmp(s, "1") == 0) : (strstr(s, "1.0") != NULL), "string (one)");
    flint_free(s);

    gr_heap_clear(x, ctx); gr_heap_clear(y, ctx); gr_heap_clear(z, ctx);
    acb_clear(a); arb_clear(t); fmpz_clear(f);
    gr_ctx_clear(ctx); gr_ctx_clear(R); gr_ctx_clear(RR); gr_ctx_clear(CC);
}

/* the typed interface (a few calls per type, as a check of the
   instantiations) */
static void
_t_c_typed(void)
{
    d2c_t x, y, z;
    d2cb_t a, b, c;
    d2_t r;
    d2b_t rb;
    int n = 2, flags = 0;

    d2c_set_d_d(x, 3.0, 4.0);
    d2c_abs(r, x);
    C_CHECK(r->d[0] == 5.0 && r->d[1] == 0.0, "d2c_abs");
    d2c_conj(y, x);
    d2c_mul(z, x, y);
    C_CHECK(z->re.d[0] == 25.0 && z->im.d[0] == 0.0, "d2c_mul");
    d2c_sqr(z, x);
    C_CHECK(z->re.d[0] == -7.0 && z->im.d[0] == 24.0, "d2c_sqr");
    d2c_sqrt(y, z);
    C_CHECK(d2c_equal(y, x), "d2c_sqrt");
    d2c_mul_onei(y, x);
    C_CHECK(y->re.d[0] == -4.0 && y->im.d[0] == 3.0, "d2c_mul_onei");

    d2cb_set_d_d(a, 3.0, 4.0);
    d2cb_abs(rb, a);
    C_CHECK(rb->d[0] == 5.0 && rb->rad == 0.0, "d2cb_abs");
    d2cb_sqr(c, a);
    C_CHECK(c->re.d[0] == -7.0 && c->im.d[0] == 24.0 && d2cb_is_exact(c), "d2cb_sqr");
    d2cb_sqrt(b, c);
    C_CHECK(d2cb_equal(b, a), "d2cb_sqrt");
    C_CHECK(d2cb_div(c, c, a) == GR_SUCCESS && d2cb_contains(c, a), "d2cb_div");
    d2cb_zero(b);
    C_CHECK(d2cb_div(c, a, b) == GR_DOMAIN && d2cb_log(c, b) == GR_DOMAIN, "d2cb zero");
}

/* the whole typed interface of dNc and dNcb, for every N, against the
   generic methods of the corresponding contexts (which call the same
   kernels: the results are identical where the method succeeds) */
#define T_SAME(z, w) (memcmp((z), (w), sizeof((z)[0])) == 0)
#define T_C1(N, f, g) \
    do { \
        d##N##c_##f(z, x); \
        if (g(w, x, C) == GR_SUCCESS) \
            C_CHECK(T_SAME(z, w), "d" #N "c_" #f); \
    } while (0)
#define T_C2(N, f, g) \
    do { \
        d##N##c_##f(z, x, y); \
        if (g(w, x, y, C) == GR_SUCCESS) \
            C_CHECK(T_SAME(z, w), "d" #N "c_" #f); \
    } while (0)
#define T_CB1(N, f, g) \
    do { \
        int _st = GR_SUCCESS; \
        _st = (int) d##N##cb_##f(c, a) + 0; \
        if (g(e, a, CB) == GR_SUCCESS) \
            C_CHECK(T_SAME(c, e), "d" #N "cb_" #f); \
        (void) _st; \
    } while (0)
#define T_CB1V(N, f, g) \
    do { \
        d##N##cb_##f(c, a); \
        if (g(e, a, CB) == GR_SUCCESS) \
            C_CHECK(T_SAME(c, e), "d" #N "cb_" #f); \
    } while (0)
#define T_CB2(N, f, g) \
    do { \
        int _st = d##N##cb_##f(c, a, b); \
        if (g(e, a, b, CB) == GR_SUCCESS) \
            C_CHECK(_st == GR_SUCCESS && T_SAME(c, e), "d" #N "cb_" #f); \
    } while (0)
#define T_CB2V(N, f, g) \
    do { \
        d##N##cb_##f(c, a, b); \
        if (g(e, a, b, CB) == GR_SUCCESS) \
            C_CHECK(T_SAME(c, e), "d" #N "cb_" #f); \
    } while (0)

#define T_TYPED_ALL(N) \
static void \
_t_c_typed_all_##N(flint_rand_t state, slong iters) \
{ \
    gr_ctx_t C, CB; \
    d##N##c_t x, y, z, w; \
    d##N##cb_t a, b, c, e; \
    d##N##_t r; \
    d##N##b_t rb; \
    acb_t t, u; \
    char * str; \
    slong it; \
    int n = N, flags = 0; \
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(C, N, DFLOAT_COMPLEX)); \
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(CB, N, DFLOAT_COMPLEX | DFLOAT_BALL)); \
    acb_init(t); \
    acb_init(u); \
    for (it = 0; it < iters; it++) \
    { \
        if (it % 4 == 0) { d##N##c_randtest_special(x, state); d##N##c_randtest(y, state); } \
        else { d##N##c_randtest(x, state); d##N##c_randtest(y, state); } \
        if (it % 4 == 1) { d##N##cb_randtest_special(a, state); d##N##cb_randtest(b, state); } \
        else { d##N##cb_randtest(a, state); d##N##cb_randtest(b, state); } \
        /* plain */ \
        T_C1(N, neg, gr_neg); T_C1(N, conj, gr_conj); T_C1(N, sqr, gr_sqr); \
        T_C1(N, inv, gr_inv); T_C1(N, sqrt, gr_sqrt); T_C1(N, rsqrt, gr_rsqrt); \
        T_C1(N, exp, gr_exp); T_C1(N, log, gr_log); T_C1(N, sin, gr_sin); \
        T_C1(N, cos, gr_cos); T_C1(N, sinh, gr_sinh); T_C1(N, cosh, gr_cosh); \
        T_C1(N, tan, gr_tan); T_C1(N, tanh, gr_tanh); \
        T_C2(N, add, gr_add); T_C2(N, sub, gr_sub); T_C2(N, mul, gr_mul); \
        T_C2(N, div, gr_div); T_C2(N, pow, gr_pow); \
        d##N##c_mul_2exp_si(z, x, 3); \
        GR_MUST_SUCCEED(gr_mul_2exp_si(w, x, 3, C)); \
        C_CHECK(T_SAME(z, w), "d" #N "c_mul_2exp_si"); \
        d##N##c_onei(w); d##N##c_mul_onei(z, x); \
        if (gr_mul(w, x, w, C) == GR_SUCCESS && d##N##c_is_finite(x)) \
            C_CHECK(d##N##c_equal(z, w), "d" #N "c_mul_onei"); \
        d##N##c_set_re_im(w, &y->re, &y->im); C_CHECK(T_SAME(w, y), "d" #N "c_set_re_im"); \
        d##N##c_mul_re(z, x, &y->re); \
        d##N##c_set_re_im(w, &y->re, &y->re); d##N##_zero(&w->im); \
        if (gr_mul(w, x, w, C) == GR_SUCCESS) \
            C_CHECK(T_SAME(z, w), "d" #N "c_mul_re"); \
        d##N##c_sin_cos(z, w, x); \
        { d##N##c_t s1, c1; d##N##c_sin(s1, x); d##N##c_cos(c1, x); \
          C_CHECK(T_SAME(z, s1) && T_SAME(w, c1), "d" #N "c_sin_cos"); } \
        d##N##c_abs(r, x); \
        if (gr_abs(w, x, C) == GR_SUCCESS) C_CHECK(T_SAME(r, &w->re), "d" #N "c_abs"); \
        d##N##c_arg(r, x); \
        if (gr_arg(w, x, C) == GR_SUCCESS) C_CHECK(T_SAME(r, &w->re), "d" #N "c_arg"); \
        d##N##c_set_d_d(z, 1.5, -0.25); d##N##c_get_acb(t, z); \
        acb_set_d_d(u, 1.5, -0.25); C_CHECK(acb_equal(t, u), "d" #N "c_set_d_d"); \
        d##N##c_set_d(z, 1.0); d##N##c_one(w); \
        C_CHECK(d##N##c_equal(z, w) && d##N##c_is_one(w) && d##N##c_is_real(w) && !d##N##c_is_zero(w), "d" #N "c_one"); \
        d##N##c_zero(z); C_CHECK(d##N##c_is_zero(z), "d" #N "c_zero"); \
        d##N##c_set(z, x); C_CHECK(T_SAME(z, x), "d" #N "c_set"); \
        if (d##N##c_is_finite(x)) \
        { \
            d##N##c_get_acb(t, x); d##N##c_set_acb(z, t); d##N##c_get_acb(u, z); \
            C_CHECK(acb_equal(t, u), "d" #N "c_set_acb"); \
        } \
        str = d##N##c_get_str(x, 10); C_CHECK(str != NULL, "d" #N "c_get_str"); flint_free(str); \
        /* balls */ \
        T_CB1V(N, neg, gr_neg); T_CB1V(N, conj, gr_conj); T_CB1V(N, sqr, gr_sqr); \
        T_CB1(N, inv, gr_inv); T_CB1V(N, sqrt, gr_sqrt); T_CB1(N, rsqrt, gr_rsqrt); \
        T_CB1V(N, exp, gr_exp); T_CB1(N, log, gr_log); T_CB1V(N, sin, gr_sin); \
        T_CB1V(N, cos, gr_cos); T_CB1V(N, sinh, gr_sinh); T_CB1V(N, cosh, gr_cosh); \
        T_CB1(N, tan, gr_tan); T_CB1(N, tanh, gr_tanh); \
        T_CB2V(N, add, gr_add); T_CB2V(N, sub, gr_sub); T_CB2V(N, mul, gr_mul); \
        T_CB2(N, div, gr_div); T_CB2(N, pow, gr_pow); \
        d##N##cb_mul_2exp_si(c, a, -2); \
        GR_MUST_SUCCEED(gr_mul_2exp_si(e, a, -2, CB)); \
        C_CHECK(T_SAME(c, e), "d" #N "cb_mul_2exp_si"); \
        d##N##cb_onei(e); d##N##cb_mul_onei(c, a); \
        if (gr_mul(e, a, e, CB) == GR_SUCCESS && d##N##cb_is_finite(a)) \
            C_CHECK(d##N##cb_overlaps(c, e), "d" #N "cb_mul_onei"); \
        d##N##cb_set_re_im(e, &b->re, &b->im); C_CHECK(T_SAME(e, b), "d" #N "cb_set_re_im"); \
        d##N##cb_mul_re(c, a, &b->re); \
        d##N##cb_set_re_im(e, &b->re, &b->re); d##N##b_zero(&e->im); \
        if (gr_mul(e, a, e, CB) == GR_SUCCESS && d##N##cb_is_finite(a) && d##N##cb_is_finite(b)) \
            C_CHECK(d##N##cb_overlaps(c, e), "d" #N "cb_mul_re"); \
        d##N##cb_sin_cos(c, e, a); \
        { d##N##cb_t s1, c1; d##N##cb_sin(s1, a); d##N##cb_cos(c1, a); \
          C_CHECK(T_SAME(c, s1) && T_SAME(e, c1), "d" #N "cb_sin_cos"); } \
        d##N##cb_abs(rb, a); \
        if (gr_abs(e, a, CB) == GR_SUCCESS) C_CHECK(T_SAME(rb, &e->re), "d" #N "cb_abs"); \
        d##N##cb_arg(rb, a); \
        if (gr_arg(e, a, CB) == GR_SUCCESS) C_CHECK(T_SAME(rb, &e->re), "d" #N "cb_arg"); \
        d##N##cb_set_d_d(c, 1.5, -0.25); \
        C_CHECK(d##N##cb_is_exact(c) && !d##N##cb_is_real(c) && !d##N##cb_contains_zero(c), "d" #N "cb_set_d_d"); \
        d##N##cb_set_d(c, 1.0); d##N##cb_one(e); \
        C_CHECK(d##N##cb_equal(c, e) && d##N##cb_is_one(e) && d##N##cb_contains(e, c), "d" #N "cb_one"); \
        d##N##cb_zero(c); C_CHECK(d##N##cb_is_zero(c) && d##N##cb_contains_zero(c), "d" #N "cb_zero"); \
        d##N##cb_indeterminate(c); C_CHECK(!d##N##cb_is_finite(c) && d##N##cb_contains(c, a), "d" #N "cb_indeterminate"); \
        d##N##cb_set(c, a); C_CHECK(T_SAME(c, a), "d" #N "cb_set"); \
        d##N##cb_get_acb(t, a); d##N##cb_set_acb(c, t); d##N##cb_get_acb(u, c); \
        C_CHECK(acb_contains(u, t), "d" #N "cb_set_acb"); \
        str = d##N##cb_get_str(a, 10); C_CHECK(str != NULL, "d" #N "cb_get_str"); flint_free(str); \
    } \
    acb_clear(t); \
    acb_clear(u); \
    gr_ctx_clear(C); \
    gr_ctx_clear(CB); \
}

T_TYPED_ALL(1)
T_TYPED_ALL(2)
T_TYPED_ALL(3)
T_TYPED_ALL(4)

#undef T_TYPED_ALL
#undef T_SAME
#undef T_C1
#undef T_C2
#undef T_CB1
#undef T_CB1V
#undef T_CB2
#undef T_CB2V

TEST_FUNCTION_START(complex, state)
{
    gr_ctx_t ctx;
    int n, flags;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    _t_c_typed();
    _t_c_typed_all_1(state, 50 * flint_test_multiplier());
    _t_c_typed_all_2(state, 50 * flint_test_multiplier());
    _t_c_typed_all_3(state, 50 * flint_test_multiplier());
    _t_c_typed_all_4(state, 50 * flint_test_multiplier());

    for (n = 1; n <= DFLOAT_MAX_N; n++)
    {
        for (flags = DFLOAT_COMPLEX; flags <= (DFLOAT_COMPLEX | 7); flags++)
        {
            int ball = (flags & DFLOAT_BALL) != 0;
            GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, flags));
            if (ball)
                gr_test_ring(ctx, 100, 0);
            else
                gr_test_floating_point(ctx, 100, 0);
            _t_c_functions(ctx, n, flags, state, 300 * flint_test_multiplier());
            _t_c_vectors(ctx, n, flags, state, 20 * flint_test_multiplier());
            gr_ctx_clear(ctx);
            _t_c_special(n, flags);
        }
    }

    TEST_FUNCTION_END(state);
}
#undef C_CHECK
