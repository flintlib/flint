/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"
#include "ulong_extras.h"

/* reference: X[k] = sum_j v[j] prod_a w_a^(j_a k_a) with w_a the
   canonical root of order cyc[a], directly from the definition */
static int
_prod_dft_naive(gr_ptr res, gr_srcptr vec, const ulong * cyc, slong num,
        gr_srcptr w, ulong order, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong n = 1, j, k;
    slong a;
    gr_ptr roots, t;

    for (a = 0; a < num; a++)
        n *= cyc[a];

    /* roots[a] = canonical root of order cyc[a] */
    roots = gr_heap_init_vec(FLINT_MAX(num, 1), ctx);
    for (a = 0; a < num; a++)
    {
        if (w == NULL)
            status |= gr_dft_default_root(GR_ENTRY(roots, a, sz),
                    cyc[a], ctx);
        else
            status |= gr_pow_ui(GR_ENTRY(roots, a, sz), w,
                    order / cyc[a], ctx);
    }

    GR_TMP_INIT(t, ctx);

    for (k = 0; k < n; k++)
    {
        gr_ptr acc = GR_ENTRY(res, k, sz);
        status |= gr_zero(acc, ctx);

        for (j = 0; j < n; j++)
        {
            ulong jj = j, kk = k, inner = n;

            status |= gr_set(t, GR_ENTRY(vec, j, sz), ctx);
            for (a = 0; a < num; a++)
            {
                ulong m = cyc[a], ja, ka;
                gr_ptr p;

                inner /= m;
                ja = (jj / inner) % m;
                ka = (kk / inner) % m;

                GR_TMP_INIT(p, ctx);
                status |= gr_pow_ui(p, GR_ENTRY(roots, a, sz),
                        (ja * ka) % m, ctx);
                status |= gr_mul(t, t, p, ctx);
                GR_TMP_CLEAR(p, ctx);
            }
            status |= gr_add(acc, acc, t, ctx);
        }
    }

    GR_TMP_CLEAR(t, ctx);
    gr_heap_clear_vec(roots, FLINT_MAX(num, 1), ctx);
    return status;
}

TEST_FUNCTION_START(gr_dft_prod, state)
{
    slong iter;

    /* deterministic sweep over ring, shape, and threading */
    {
        static const ulong shapes[][5] = {
            { 2, 0 }, { 4, 0 }, { 2, 3, 0 }, { 3, 5, 0 }, { 2, 2, 2, 2, 0 },
            { 4, 4, 0 }, { 6, 5, 0 }, { 1, 3, 1, 0 }, { 7, 0 }, { 8, 3, 0 },
        };
        slong si, ring, threads;

        for (ring = 0; ring < 2; ring++)
        for (si = 0; si < (slong) (sizeof(shapes) / sizeof(shapes[0])); si++)
        for (threads = 1; threads <= 3; threads += 2)
        {
            gr_ctx_t ctx;
            ulong cyc[5], n = 1, j;
            slong num = 0, sz;
            gr_ptr x, y, z;
            gr_dft_prod_pre_t P;
            int status = GR_SUCCESS;

            while (shapes[si][num] != 0)
            {
                cyc[num] = shapes[si][num];
                n *= cyc[num];
                num++;
            }

            flint_set_num_threads(threads);

            if (ring == 0)
                gr_ctx_init_complex_acb(ctx, 128);
            else
            {
                /* prime q = 1 mod lcm: q = 1 mod 840 covers all shapes */
                gr_ctx_init_nmod(ctx, 68041);   /* 68041 = 81 * 840 + 1 */
            }
            sz = ctx->sizeof_elem;

            x = gr_heap_init_vec(n, ctx);
            y = gr_heap_init_vec(n, ctx);
            z = gr_heap_init_vec(n, ctx);
            for (j = 0; j < n; j++)
                status |= gr_set_ui(GR_ENTRY(x, j, sz),
                        1 + (j * j * 17 + j) % 97, ctx);

            if (ring == 0)
            {
                status |= _prod_dft_naive(y, x, cyc, num, NULL, 0, ctx);
                status |= gr_dft_prod_precomp_init(P, cyc, num, 0, ctx);
            }
            else
            {
                /* w = generator^((q-1)/840): a root of order 840,
                   divisible by every component length in the sweep */
                gr_ptr w;
                GR_TMP_INIT(w, ctx);
                status |= gr_set_ui(w, n_primitive_root_prime(68041), ctx);
                status |= gr_pow_ui(w, w, (68041 - 1) / 840, ctx);
                status |= _prod_dft_naive(y, x, cyc, num, w, 840, ctx);
                status |= gr_dft_prod_precomp_init_root(P, w, 840, cyc,
                        num, 0, ctx);
                GR_TMP_CLEAR(w, ctx);
            }
            gr_dft_prod_precomp_set_serial_block(P, 1);
            status |= gr_dft_prod_precomp(z, x, P, ctx);

            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("status, ring %wd shape %wd\n", ring, si);

            for (j = 0; j < n; j++)
                if (gr_equal(GR_ENTRY(z, j, sz), GR_ENTRY(y, j, sz), ctx)
                        == T_FALSE)
                    TEST_FUNCTION_FAIL("forward, ring %wd shape %wd j %wu\n",
                            ring, si, j);

            /* aliased forward, then inverse roundtrip */
            status |= _gr_vec_set(z, x, n, ctx);
            status |= gr_dft_prod_precomp(z, z, P, ctx);
            status |= gr_dft_prod_inverse_precomp(z, z, P, ctx);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("roundtrip status, shape %wd\n", si);
            for (j = 0; j < n; j++)
                if (gr_equal(GR_ENTRY(z, j, sz), GR_ENTRY(x, j, sz), ctx)
                        == T_FALSE)
                    TEST_FUNCTION_FAIL("roundtrip, ring %wd shape %wd j %wu\n",
                            ring, si, j);

            gr_dft_prod_precomp_clear(P);
            gr_heap_clear_vec(x, n, ctx);
            gr_heap_clear_vec(y, n, ctx);
            gr_heap_clear_vec(z, n, ctx);
            gr_ctx_clear(ctx);
            flint_set_num_threads(1);
        }
    }

    /* acb wrapper: fixed-point vs ball containment (both enclose the
       exact transform, so they must overlap), plus roundtrip */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        ulong cyc[4], n = 1, j;
        slong num = 1 + n_randint(state, 3), a;
        slong prec = 32 + n_randint(state, 300);
        acb_ptr v, w1, w2;

        for (a = 0; a < num; a++)
        {
            cyc[a] = 1 + n_randint(state, 8);
            n *= cyc[a];
        }

        v = _acb_vec_init(n);
        w1 = _acb_vec_init(n);
        w2 = _acb_vec_init(n);

        for (j = 0; j < n; j++)
        {
            acb_randtest(v + j, state, prec, 4);
            if (n_randint(state, 8) == 0)
                acb_mul_2exp_si(v + j, v + j,
                        (slong) n_randint(state, 200) - 100);
        }

        if (_gr_dft_acb_prod(w1, v, cyc, num, 0, 2, prec) == GR_SUCCESS)
        {
            if (_gr_dft_acb_prod(w2, v, cyc, num, 0, 1, prec) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("ball path failed\n");
            for (j = 0; j < n; j++)
                if (!acb_overlaps(w1 + j, w2 + j))
                    TEST_FUNCTION_FAIL("containment, n %wu j %wu\n", n, j);
        }

        /* roundtrip through the public entries */
        gr_dft_acb_prod(w1, v, cyc, num, prec);
        gr_dft_acb_prod_inverse(w2, w1, cyc, num, prec);
        for (j = 0; j < n; j++)
            if (!acb_contains(w2 + j, v + j) && !acb_overlaps(w2 + j, v + j))
                TEST_FUNCTION_FAIL("acb roundtrip, n %wu j %wu\n", n, j);

        /* plan interface must agree with the one-shot */
        {
            gr_dft_acb_prod_pre_t Q;
            if (gr_dft_acb_prod_precomp_init(Q, cyc, num, prec)
                    == GR_SUCCESS)
            {
                gr_dft_acb_prod_precomp(w2, v, Q, prec);
                for (j = 0; j < n; j++)
                    if (!acb_overlaps(w1 + j, w2 + j))
                        TEST_FUNCTION_FAIL("plan vs one-shot, j %wu\n", j);
                gr_dft_acb_prod_precomp_clear(Q);
            }
        }

        _acb_vec_clear(v, n);
        _acb_vec_clear(w1, n);
        _acb_vec_clear(w2, n);
    }

    TEST_FUNCTION_END(state);
}
