/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Targeted tests for special cases of the fq, fq_nmod and fq_zech
   methods: multiplication by powers of two (characteristic 2, c == 0,
   aliasing), dot products (degree above the in-place cutoff, lengths
   0 and 1, aliasing with the initial value) and fq_nmod multiplication
   by an ulong into an element with small allocation. */

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "gr.h"
#include "gr_vec.h"

/* reference: res = x * 2^c by repeated multiplication or division by 2 */
static int
_fq_special_mul_2exp_ref(gr_ptr res, gr_srcptr x, slong c, gr_ctx_t ctx)
{
    slong i;
    int status = gr_set(res, x, ctx);

    for (i = 0; i < FLINT_ABS(c); i++)
    {
        if (c > 0)
            status |= gr_mul_ui(res, res, 2, ctx);
        else
            status |= gr_div_ui(res, res, 2, ctx);
    }

    return status;
}

static int
_fq_special_dot_ref(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, int rev, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr s, t;

    s = gr_heap_init(ctx);
    t = gr_heap_init(ctx);

    if (initial != NULL)
        status |= gr_set(s, initial, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, rev ? len - 1 - i : i, sz), ctx);
        if (subtract)
            status |= gr_sub(s, s, t, ctx);
        else
            status |= gr_add(s, s, t, ctx);
    }

    gr_swap(res, s, ctx);

    gr_heap_clear(s, ctx);
    gr_heap_clear(t, ctx);

    return status;
}

static void
_fq_special_check_ctx(gr_ctx_t ctx, int char2, flint_rand_t state)
{
    slong sz = ctx->sizeof_elem;
    slong c, len, i, j;
    int status, expected, alias, subtract, rev, init_kind;
    gr_ptr x, r, s, vec, rvec, v1, v2, init;

    x = gr_heap_init(ctx);
    r = gr_heap_init(ctx);
    s = gr_heap_init(ctx);
    init = gr_heap_init(ctx);
    vec = gr_heap_init_vec(6, ctx);
    rvec = gr_heap_init_vec(6, ctx);
    v1 = gr_heap_init_vec(6, ctx);
    v2 = gr_heap_init_vec(6, ctx);

    /* mul_2exp_si and vec_mul_scalar_2exp_si */
    for (c = -3; c <= 3; c++)
    {
        expected = (char2 && c < 0) ? GR_DOMAIN : GR_SUCCESS;

        for (alias = 0; alias < 2; alias++)
        {
            GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
            GR_MUST_SUCCEED(gr_randtest(r, state, ctx));

            if (alias)
            {
                GR_MUST_SUCCEED(gr_set(r, x, ctx));
                status = gr_mul_2exp_si(r, r, c, ctx);
            }
            else
            {
                status = gr_mul_2exp_si(r, x, c, ctx);
            }

            if (status != expected)
            {
                flint_printf("FAIL: mul_2exp_si status\n");
                gr_ctx_println(ctx);
                flint_printf("c = %wd, alias = %d, status = %d\n", c, alias, status);
                flint_abort();
            }

            if (status == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(_fq_special_mul_2exp_ref(s, x, c, ctx));

                if (gr_equal(r, s, ctx) != T_TRUE ||
                    (char2 && c > 0 && gr_is_zero(r, ctx) != T_TRUE) ||
                    (c == 0 && gr_equal(r, x, ctx) != T_TRUE))
                {
                    flint_printf("FAIL: mul_2exp_si value\n");
                    gr_ctx_println(ctx);
                    flint_printf("c = %wd, alias = %d\n", c, alias);
                    flint_printf("x = %{gr}\nr = %{gr}\ns = %{gr}\n", x, ctx, r, ctx, s, ctx);
                    flint_abort();
                }
            }

            for (len = 0; len <= 5; len++)
            {
                GR_MUST_SUCCEED(_gr_vec_randtest(vec, state, len, ctx));
                GR_MUST_SUCCEED(_gr_vec_randtest(rvec, state, len, ctx));

                if (alias)
                {
                    GR_MUST_SUCCEED(_gr_vec_set(rvec, vec, len, ctx));
                    status = _gr_vec_mul_scalar_2exp_si(rvec, rvec, len, c, ctx);
                }
                else
                {
                    status = _gr_vec_mul_scalar_2exp_si(rvec, vec, len, c, ctx);
                }

                /* a zero-length vector in characteristic 2 with c < 0
                   may or may not flag an error */
                if (status != expected && !(len == 0 && status == GR_SUCCESS))
                {
                    flint_printf("FAIL: vec_mul_scalar_2exp_si status\n");
                    gr_ctx_println(ctx);
                    flint_printf("c = %wd, len = %wd, alias = %d, status = %d\n", c, len, alias, status);
                    flint_abort();
                }

                if (status == GR_SUCCESS)
                {
                    for (i = 0; i < len; i++)
                    {
                        GR_MUST_SUCCEED(_fq_special_mul_2exp_ref(s, GR_ENTRY(vec, i, sz), c, ctx));

                        if (gr_equal(GR_ENTRY(rvec, i, sz), s, ctx) != T_TRUE)
                        {
                            flint_printf("FAIL: vec_mul_scalar_2exp_si value\n");
                            gr_ctx_println(ctx);
                            flint_printf("c = %wd, len = %wd, alias = %d, i = %wd\n", c, len, alias, i);
                            flint_printf("vec = %{gr*}\nrvec = %{gr*}\n", vec, len, ctx, rvec, len, ctx);
                            flint_abort();
                        }
                    }
                }
            }
        }
    }

    /* dot products */
    for (len = 0; len <= 5; len++)
    {
        for (j = 0; j < 16; j++)
        {
            subtract = j & 1;
            rev = (j >> 1) & 1;
            init_kind = (j >> 2);   /* 0: NULL, 1: initial, 2: res == initial, 3: NULL, zero res */

            GR_MUST_SUCCEED(_gr_vec_randtest(v1, state, len, ctx));
            GR_MUST_SUCCEED(_gr_vec_randtest(v2, state, len, ctx));
            GR_MUST_SUCCEED(gr_randtest(init, state, ctx));

            /* output initialised to elements of various lengths */
            if (init_kind == 3)
                GR_MUST_SUCCEED(gr_zero(r, ctx));
            else
                GR_MUST_SUCCEED(gr_randtest(r, state, ctx));

            GR_MUST_SUCCEED(_fq_special_dot_ref(s, (init_kind == 1 || init_kind == 2) ? init : NULL,
                subtract, v1, v2, len, rev, ctx));

            if (init_kind == 2)
            {
                GR_MUST_SUCCEED(gr_set(r, init, ctx));
                if (rev)
                    status = _gr_vec_dot_rev(r, r, subtract, v1, v2, len, ctx);
                else
                    status = _gr_vec_dot(r, r, subtract, v1, v2, len, ctx);
            }
            else
            {
                if (rev)
                    status = _gr_vec_dot_rev(r, init_kind == 1 ? init : NULL, subtract, v1, v2, len, ctx);
                else
                    status = _gr_vec_dot(r, init_kind == 1 ? init : NULL, subtract, v1, v2, len, ctx);
            }

            if (status != GR_SUCCESS || gr_equal(r, s, ctx) != T_TRUE)
            {
                flint_printf("FAIL: dot\n");
                gr_ctx_println(ctx);
                flint_printf("len = %wd, subtract = %d, rev = %d, init_kind = %d, status = %d\n",
                    len, subtract, rev, init_kind, status);
                flint_printf("v1 = %{gr*}\nv2 = %{gr*}\n", v1, len, ctx, v2, len, ctx);
                flint_printf("init = %{gr}\nr = %{gr}\ns = %{gr}\n", init, ctx, r, ctx, s, ctx);
                flint_abort();
            }
        }
    }

    gr_heap_clear(x, ctx);
    gr_heap_clear(r, ctx);
    gr_heap_clear(s, ctx);
    gr_heap_clear(init, ctx);
    gr_heap_clear_vec(vec, 6, ctx);
    gr_heap_clear_vec(rvec, 6, ctx);
    gr_heap_clear_vec(v1, 6, ctx);
    gr_heap_clear_vec(v2, 6, ctx);
}

TEST_FUNCTION_START(gr_fq_special, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        ulong p;
        slong d;
        int char2, which;

        char2 = n_randint(state, 2);
        which = n_randint(state, 3);

        if (which == 2)
        {
            /* fq_zech: keep the field small */
            if (char2)
            {
                p = 2;
                d = 1 + n_randint(state, 12);
            }
            else
            {
                p = n_randprime(state, 2 + n_randint(state, 4), 1);
                if (p == 2)
                    p = 3;
                d = 1 + n_randint(state, (p < 8) ? 5 : 2);
            }

            gr_ctx_init_fq_zech(ctx, p, d, "a");
        }
        else
        {
            if (char2)
                p = 2;
            else
            {
                p = n_randtest_prime(state, 0);
                if (p == 2)
                    p = 3;
            }

            d = 1 + n_randint(state, 12);

            if (which == 0)
            {
                fmpz_t P;
                fmpz_init(P);

                if (!char2 && n_randint(state, 4) == 0)
                    fmpz_randprime(P, state, FLINT_BITS + 1 + n_randint(state, 40), 0);
                else
                    fmpz_set_ui(P, p);

                gr_ctx_init_fq(ctx, P, d, "a");
                fmpz_clear(P);
            }
            else
            {
                gr_ctx_init_fq_nmod(ctx, p, d, "a");
            }
        }

        _fq_special_check_ctx(ctx, char2, state);

        /* fq_nmod multiplication by an ulong, with an output element
           that has less space allocated than required */
        if (which == 1)
        {
          slong k;
          for (k = 0; k < 5; k++)
          {
            nmod_poly_t res;
            gr_ptr x, s;
            ulong y;

            x = gr_heap_init(ctx);
            s = gr_heap_init(ctx);

            GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
            y = n_randtest(state);

            nmod_poly_init(res, p);
            GR_MUST_SUCCEED(gr_mul_ui((gr_ptr) res, x, y, ctx));

            GR_MUST_SUCCEED(gr_set_ui(s, y, ctx));
            GR_MUST_SUCCEED(gr_mul(s, x, s, ctx));

            if (gr_equal((gr_ptr) res, s, ctx) != T_TRUE)
            {
                flint_printf("FAIL: fq_nmod mul_ui\n");
                gr_ctx_println(ctx);
                flint_printf("y = %wu\nx = %{gr}\ns = %{gr}\n", y, x, ctx, s, ctx);
                flint_abort();
            }

            nmod_poly_clear(res);
            gr_heap_clear(x, ctx);
            gr_heap_clear(s, ctx);
          }
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
