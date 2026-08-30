/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "test_helpers.h"
#include "gr_ore_poly.h"

static void
check_gen_action(gr_srcptr x, gr_srcptr expected, gr_ore_poly_ctx_t octx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_ore_poly_t D;
    gr_ptr got;

    gr_ore_poly_init(D, octx);
    got = gr_heap_init(cctx);

    status |= gr_ore_poly_gen(D, octx);
    status |= gr_ore_poly_apply(got, D, x, octx);

    if (status != GR_SUCCESS || gr_equal(got, expected, cctx) != T_TRUE)
    {
        flint_printf("FAIL: D(gen)\n\noctx = %{gr_ctx}\n", octx);
        flint_abort();
    }

    gr_heap_clear(got, cctx);
    gr_ore_poly_clear(D, octx);
}

static void
check_gen_actions(flint_rand_t state)
{
    int status = GR_SUCCESS;
    gr_ctx_t zctx, cctx;
    gr_ore_poly_ctx_t octx;
    gr_ptr x, q, expected;

    gr_ctx_init_fmpz(zctx);
    gr_ctx_init_gr_poly(cctx, zctx);
    x = gr_heap_init(cctx);
    q = gr_heap_init(cctx);
    expected = gr_heap_init(cctx);
    status |= gr_gen(x, cctx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_DERIVATIVE);
    status |= gr_one(expected, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
    status |= gr_set(expected, x, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_FORWARD_SHIFT);
    status |= gr_add_si(expected, x, 1, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_FORWARD_DIFFERENCE);
    status |= gr_one(expected, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_BACKWARD_SHIFT);
    status |= gr_add_si(expected, x, -1, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_ore_poly_ctx_init(octx, cctx, 0, ORE_ALGEBRA_BACKWARD_DIFFERENCE);
    status |= gr_one(expected, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    status |= gr_set_si(q, 3, cctx);
    status |= gr_ore_poly_ctx_init_q_shift(octx, cctx, 0, q);
    status |= gr_mul(expected, q, x, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    status |= gr_ore_poly_ctx_init_mahler(octx, cctx, 0, 3);
    status |= gr_pow_ui(expected, x, 3, cctx);
    check_gen_action(x, expected, octx, cctx);
    gr_ore_poly_ctx_clear(octx);

    gr_heap_clear(x, cctx);
    gr_heap_clear(q, cctx);
    gr_heap_clear(expected, cctx);
    gr_ctx_clear(cctx);
    gr_ctx_clear(zctx);

    {
        fmpz_t p;
        gr_ctx_t fctx;
        gr_ptr fx, fexpected;

        fmpz_init_set_ui(p, 17);
        gr_ctx_init_fq(fctx, p, 3, NULL);
        fmpz_clear(p);

        fx = gr_heap_init(fctx);
        fexpected = gr_heap_init(fctx);

        status |= gr_gen(fx, fctx);
        status |= gr_pow_ui(fexpected, fx, 17, fctx);

        gr_ore_poly_ctx_init(octx, fctx, 0, ORE_ALGEBRA_FROBENIUS);
        check_gen_action(fx, fexpected, octx, fctx);
        gr_ore_poly_ctx_clear(octx);

        gr_heap_clear(fx, fctx);
        gr_heap_clear(fexpected, fctx);
        gr_ctx_clear(fctx);
    }

    if (status != GR_SUCCESS)
    {
        flint_printf("FAIL: unexpected failure\n");
        flint_abort();
    }
}

TEST_GR_FUNCTION_START(gr_ore_poly_apply, state, count_success, count_domain, count_unable)
{
    check_gen_actions(state);

    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        slong maxlen;
        int status = GR_SUCCESS;

        gr_ore_poly_ctx_init_randtest2(cctx, ctx, state);

        if (GR_ORE_POLY_CTX(ctx)->which_algebra == ORE_ALGEBRA_MAHLER
            || GR_ORE_POLY_CTX(ctx)->which_algebra == ORE_ALGEBRA_Q_SHIFT)
            maxlen = 2;
        else
            maxlen = 5;

        gr_ore_poly_t P, Q, PQ, PpQ, D;
        gr_ore_poly_init(P, ctx);
        gr_ore_poly_init(Q, ctx);
        gr_ore_poly_init(PQ, ctx);
        gr_ore_poly_init(PpQ, ctx);
        gr_ore_poly_init(D, ctx);

        gr_ptr f = gr_heap_init(cctx);
        gr_ptr g = gr_heap_init(cctx);
        gr_ptr d1 = gr_heap_init(cctx);
        gr_ptr u = gr_heap_init(cctx);
        gr_ptr v = gr_heap_init(cctx);
        gr_ptr w = gr_heap_init(cctx);
        gr_ptr lhs = gr_heap_init(cctx);
        gr_ptr rhs = gr_heap_init(cctx);
        gr_ptr sf = gr_heap_init(cctx);
        gr_ptr df = gr_heap_init(cctx);
        gr_ptr one = gr_heap_init(cctx);
        gr_ptr c = gr_heap_init(cctx);

        status |= gr_ore_poly_randtest(P, state, 1 + n_randint(state, maxlen), ctx);
        status |= gr_ore_poly_randtest(Q, state, 1 + n_randint(state, maxlen), ctx);
        status |= gr_ore_poly_gen(D, ctx);
        if (n_randint(state, 8))
            status |= gr_randtest_not_zero(f, state, cctx);
        else
            status |= gr_randtest(f, state, cctx);
        status |= gr_randtest(g, state, cctx);

        status |= gr_one(one, cctx);
        status |= gr_ore_poly_apply(c, D, one, ctx);
        status |= gr_ore_poly_apply(lhs, P, f, ctx);
        status |= gr_ore_poly_apply_custom(rhs, P, f, c, ctx);
        if (status == GR_SUCCESS && gr_equal(lhs, rhs, cctx) == T_FALSE)
        {
            flint_printf("FAIL: apply(P, f) = apply_custom(P, f, D(1))\n");
            flint_abort();
        }

        if (status == GR_SUCCESS && n_randint(state, 2))
            status |= gr_set(d1, c, cctx);
        else
            status |= gr_randtest(d1, state, cctx);

        status |= gr_ore_poly_apply_custom(u, D, f, d1, ctx);
        status |= gr_ore_poly_sigma_delta(sf, df, f, ctx);
        status |= gr_mul(sf, sf, d1, cctx);
        status |= gr_add(v, sf, df, cctx);
        if (status == GR_SUCCESS && gr_equal(u, v, cctx) == T_FALSE)
        {
            flint_printf("FAIL: D(f) = sigma(f)*d1 + delta(f)\n");
            flint_abort();
        }

        status |= gr_ore_poly_mul(PQ, P, Q, ctx);
        status |= gr_ore_poly_apply_custom(lhs, PQ, f, d1, ctx);
        status |= gr_ore_poly_apply_custom(u, Q, f, d1, ctx);
        status |= gr_ore_poly_apply_custom(rhs, P, u, d1, ctx);
        if (status == GR_SUCCESS && gr_equal(lhs, rhs, cctx) == T_FALSE)
        {
            flint_printf("FAIL: (P*Q)(f) = P(Q(f))\n");
            flint_abort();
        }

        status |= gr_ore_poly_add(PpQ, P, Q, ctx);
        status |= gr_ore_poly_apply_custom(lhs, PpQ, f, d1, ctx);
        status |= gr_ore_poly_apply_custom(u, P, f, d1, ctx);
        status |= gr_ore_poly_apply_custom(v, Q, f, d1, ctx);
        status |= gr_add(rhs, u, v, cctx);
        if (status == GR_SUCCESS && gr_equal(lhs, rhs, cctx) == T_FALSE)
        {
            flint_printf("FAIL: (P+Q)(f) = P(f) + Q(f)\n");
            flint_abort();
        }

        status |= gr_add(w, f, g, cctx);
        status |= gr_ore_poly_apply_custom(lhs, P, w, d1, ctx);
        status |= gr_ore_poly_apply_custom(u, P, f, d1, ctx);
        status |= gr_ore_poly_apply_custom(v, P, g, d1, ctx);
        status |= gr_add(rhs, u, v, cctx);
        if (status == GR_SUCCESS && gr_equal(lhs, rhs, cctx) == T_FALSE)
        {
            flint_printf("FAIL: P(f+g) = P(f) + P(g)\n");
            flint_abort();
        }

        status |= gr_ore_poly_apply_custom(u, P, f, d1, ctx);
        status |= gr_set(v, f, cctx);
        status |= gr_ore_poly_apply_custom(v, P, v, d1, ctx);
        if (status == GR_SUCCESS && gr_equal(u, v, cctx) == T_FALSE)
        {
            flint_printf("FAIL: aliasing res == f\n");
            flint_abort();
        }

        status |= gr_ore_poly_apply_custom(u, P, f, d1, ctx);
        status |= gr_set(v, d1, cctx);
        status |= gr_ore_poly_apply_custom(v, P, f, v, ctx);
        if (status == GR_SUCCESS && gr_equal(u, v, cctx) == T_FALSE)
        {
            flint_printf("FAIL: aliasing res == d1\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_heap_clear(f, cctx);
        gr_heap_clear(g, cctx);
        gr_heap_clear(d1, cctx);
        gr_heap_clear(u, cctx);
        gr_heap_clear(v, cctx);
        gr_heap_clear(w, cctx);
        gr_heap_clear(lhs, cctx);
        gr_heap_clear(rhs, cctx);
        gr_heap_clear(sf, cctx);
        gr_heap_clear(df, cctx);
        gr_heap_clear(one, cctx);
        gr_heap_clear(c, cctx);
        gr_ore_poly_clear(P, ctx);
        gr_ore_poly_clear(Q, ctx);
        gr_ore_poly_clear(PQ, ctx);
        gr_ore_poly_clear(PpQ, ctx);
        gr_ore_poly_clear(D, ctx);
        gr_ore_poly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
