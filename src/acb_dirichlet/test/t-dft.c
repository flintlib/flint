/*
    Copyright (C) 2016 Pascal Molin
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
#include "acb_dirichlet.h"

/* acb_dirichlet_dft and acb_dirichlet_dft_index are thin wrappers
   around gr_dft, whose test suite validates the transforms against
   independent references; here we verify that the wrappers behave as
   documented: dft_index computes the pairing sum on a small group,
   and dft is its conjugation by the Conrey gather/scatter. */

TEST_FUNCTION_START(acb_dirichlet_dft, state)
{
    slong k;
    ulong q[6] = { 1, 2, 5, 8, 12, 30 };

    for (k = 0; k < 6; k++)
    {
        slong prec = 100;
        dirichlet_group_t G;
        dirichlet_char_t x, y;
        acb_ptr v, vq, w1, w2, t1;
        acb_t chiy, u;
        ulong i, j, len;

        dirichlet_group_init(G, q[k]);
        len = G->phi_q;

        v = _acb_vec_init(len);
        vq = _acb_vec_init(q[k]);
        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(FLINT_MAX(q[k], len));
        t1 = _acb_vec_init(len);
        acb_init(chiy);
        acb_init(u);
        dirichlet_char_init(x, G);
        dirichlet_char_init(y, G);

        for (i = 0; i < len; i++)
            acb_randtest_precise(v + i, state, prec, 0);
        for (i = 0; i < q[k]; i++)
            acb_randtest_precise(vq + i, state, prec, 0);

        /* dft_index against the pairing sum
           w[x] = sum_y conj(chi_x(y)) v[y] */
        dirichlet_char_one(x, G);
        for (i = 0; i < len; i++)
        {
            acb_zero(w1 + i);
            dirichlet_char_one(y, G);
            for (j = 0; j < len; j++)
            {
                arb_set_si(acb_realref(u),
                        -2 * (slong) dirichlet_pairing_char(G, x, y));
                arb_div_ui(acb_realref(u), acb_realref(u), G->expo,
                        prec + 10);
                arb_zero(acb_imagref(u));
                acb_exp_pi_i(chiy, u, prec + 10);
                acb_addmul(w1 + i, chiy, v + j, prec);
                dirichlet_char_next(y, G);
            }
            dirichlet_char_next(x, G);
        }

        acb_dirichlet_dft_index(w2, v, G, prec);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w1 + i, w2 + i))
                TEST_FUNCTION_FAIL("dft_index, q = %wu i = %wu\n", q[k], i);

        /* dft is dft_index conjugated by the Conrey gather/scatter */
        if (q[k] > 1)
        {
            dirichlet_char_one(x, G);
            for (i = 0; i < len; i++)
            {
                acb_set(t1 + i, vq + x->n);
                dirichlet_char_next(x, G);
            }
            acb_dirichlet_dft_index(w1, t1, G, prec);

            acb_dirichlet_dft(w2, vq, G, prec);
            dirichlet_char_one(x, G);
            for (i = 0; i < len; i++)
            {
                if (!acb_overlaps(w1 + i, w2 + x->n))
                    TEST_FUNCTION_FAIL("dft, q = %wu i = %wu\n", q[k], i);
                dirichlet_char_next(x, G);
            }
        }
        else
        {
            /* q = 1: identity on the single entry */
            acb_dirichlet_dft(w2, vq, G, prec);
            if (!acb_equal(w2, vq))
                TEST_FUNCTION_FAIL("dft, q = 1\n");
        }

        _acb_vec_clear(v, len);
        _acb_vec_clear(vq, q[k]);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, FLINT_MAX(q[k], len));
        _acb_vec_clear(t1, len);
        acb_clear(chiy);
        acb_clear(u);
        dirichlet_char_clear(x);
        dirichlet_char_clear(y);
        dirichlet_group_clear(G);
    }

    TEST_FUNCTION_END(state);
}
