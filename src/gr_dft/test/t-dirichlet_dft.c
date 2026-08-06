/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "dirichlet.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Independent naive reference (adapted from acb_dirichlet's t-dft.c,
   but not using acb_dft or acb_dirichlet, which are to become
   wrappers around this module):

       w[x] = sum_y conj(chi_x(y)) v[y]

   over the characters in lexicographic Conrey order, with the
   pairing values chi_x(y) = e(dirichlet_pairing_char(G, x, y) /
   G->expo) evaluated through a directly computed table of roots of
   unity. */
static void
_dirichlet_dft_index_naive(acb_ptr w, acb_srcptr v,
        const dirichlet_group_t G, slong prec)
{
    slong i, j, len = G->phi_q;
    ulong expo = G->expo, k;
    acb_ptr roots;
    acb_t t;
    dirichlet_char_t x, y;

    /* roots[k] = e(-k / expo), the conjugate pairing values */
    roots = _acb_vec_init(expo);
    acb_init(t);
    for (k = 0; k < expo; k++)
    {
        arb_set_si(acb_realref(t), -2 * (slong) k);
        arb_div_ui(acb_realref(t), acb_realref(t), expo, prec + 10);
        arb_zero(acb_imagref(t));
        acb_exp_pi_i(roots + k, t, prec + 10);
    }

    dirichlet_char_init(x, G);
    dirichlet_char_init(y, G);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        acb_zero(w + i);
        dirichlet_char_one(y, G);
        for (j = 0; j < len; j++)
        {
            acb_addmul(w + i,
                    roots + dirichlet_pairing_char(G, x, y), v + j, prec);
            dirichlet_char_next(y, G);
        }
        dirichlet_char_next(x, G);
    }

    dirichlet_char_clear(x);
    dirichlet_char_clear(y);
    _acb_vec_clear(roots, expo);
    acb_clear(t);
}

TEST_FUNCTION_START(gr_dft_dirichlet, state)
{
    /* moduli covering odd, even (q_even = 2, 4 and >= 8, the last
       giving two even Conrey components), prime, prime power and
       composite structure */
    slong k;
    slong nq = 16;
    ulong q[16] = { 2, 3, 4, 5, 6, 8, 10, 15, 16, 23, 24, 30, 59, 308,
        335, 961 };

    for (k = 0; k < nq; k++)
    {
        slong prec = 100;
        dirichlet_group_t G;
        acb_ptr v, vq, w1, w2;
        ulong i, len;
        gr_ctx_t ctx;
        dirichlet_char_t x;
        int status = GR_SUCCESS;

        dirichlet_group_init(G, q[k]);
        len = G->phi_q;

        v = _acb_vec_init(len);
        vq = _acb_vec_init(q[k]);
        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(FLINT_MAX(q[k], len));

        for (i = 0; i < len; i++)
            acb_randtest_precise(v + i, state, prec, 0);
        for (i = 0; i < q[k]; i++)
            acb_randtest_precise(vq + i, state, prec, 0);

        _dirichlet_dft_index_naive(w1, v, G, prec);

        /* conrey-indexed transform: the generic gr entry point on an
           acb context, and the fixed-point wrapper */
        gr_ctx_init_complex_acb(ctx, prec);
        status |= gr_dft_dirichlet_index((gr_ptr) w2, (gr_srcptr) v, G, ctx);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("index status, q = %wu\n", q[k]);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w1 + i, w2 + i))
                TEST_FUNCTION_FAIL("generic index vs naive, q = %wu "
                        "i = %wu\n", q[k], i);
        gr_ctx_clear(ctx);

        gr_dft_acb_dirichlet_index(w2, v, G, prec);
        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(w1 + i, w2 + i))
                TEST_FUNCTION_FAIL("acb index vs naive, q = %wu i = %wu\n",
                        q[k], i);
            if (acb_rel_accuracy_bits(w2 + i) < 30 && len > 1)
                TEST_FUNCTION_FAIL("index accuracy, q = %wu i = %wu "
                        "(%wd bits)\n", q[k], i,
                        acb_rel_accuracy_bits(w2 + i));
        }

        /* number-indexed transform against its definition: gather the
           values in Conrey order, transform in index order, scatter */
        {
            acb_ptr t1;

            t1 = _acb_vec_init(len);
            dirichlet_char_init(x, G);

            dirichlet_char_one(x, G);
            for (i = 0; i < len; i++)
            {
                acb_set(t1 + i, vq + x->n);
                dirichlet_char_next(x, G);
            }
            _dirichlet_dft_index_naive(w1, t1, G, prec);

            gr_ctx_init_complex_acb(ctx, prec);
            status |= gr_dft_dirichlet((gr_ptr) w2, (gr_srcptr) vq, G, ctx);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("status, q = %wu\n", q[k]);
            gr_ctx_clear(ctx);

            dirichlet_char_one(x, G);
            for (i = 0; i < len; i++)
            {
                if (!acb_overlaps(w1 + i, w2 + x->n))
                    TEST_FUNCTION_FAIL("generic number-indexed, q = %wu "
                            "i = %wu\n", q[k], i);
                dirichlet_char_next(x, G);
            }

            gr_dft_acb_dirichlet(w2, vq, G, prec);
            dirichlet_char_one(x, G);
            for (i = 0; i < len; i++)
            {
                if (!acb_overlaps(w1 + i, w2 + x->n))
                    TEST_FUNCTION_FAIL("acb number-indexed, q = %wu "
                            "i = %wu\n", q[k], i);
                dirichlet_char_next(x, G);
            }

            dirichlet_char_clear(x);
            _acb_vec_clear(t1, len);
        }

        _acb_vec_clear(v, len);
        _acb_vec_clear(vq, q[k]);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, FLINT_MAX(q[k], len));
        dirichlet_group_clear(G);
    }

    /* q = 1: the transforms reduce to the identity on the single
       entry (the Conrey number of the trivial character is 1, so the
       generic gather/scatter must not run) */
    {
        dirichlet_group_t G;
        acb_t a, b;
        gr_ctx_t ctx;
        int status = GR_SUCCESS;

        dirichlet_group_init(G, 1);
        acb_init(a);
        acb_init(b);
        acb_set_d_d(a, 1.25, -0.5);

        gr_ctx_init_complex_acb(ctx, 64);
        status |= gr_dft_dirichlet((gr_ptr) b, (gr_srcptr) a, G, ctx);
        if (status != GR_SUCCESS || !acb_equal(a, b))
            TEST_FUNCTION_FAIL("q = 1 generic\n");
        status |= gr_dft_dirichlet_index((gr_ptr) b, (gr_srcptr) a, G, ctx);
        if (status != GR_SUCCESS || !acb_equal(a, b))
            TEST_FUNCTION_FAIL("q = 1 generic index\n");
        gr_ctx_clear(ctx);

        acb_zero(b);
        gr_dft_acb_dirichlet(b, a, G, 64);
        if (!acb_equal(a, b))
            TEST_FUNCTION_FAIL("q = 1 acb\n");

        acb_zero(b);
        gr_dft_acb_dirichlet_index(b, a, G, 64);
        if (!acb_equal(a, b))
            TEST_FUNCTION_FAIL("q = 1 index\n");

        acb_clear(a);
        acb_clear(b);
        dirichlet_group_clear(G);
    }

    TEST_FUNCTION_END(state);
}
