/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_fmpq_afe_jet, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t x;
        acb_t s;
        acb_ptr t, u;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, k;
        slong prec, len, j;
        int deflate;
        mag_t mx, r;

        q = 1 + n_randint(state, n_randint(state, 4) ? 30 : 200);
        prec = 2 + n_randint(state, n_randint(state, 4) ? 200 : 600);
        len = 1 + n_randint(state, n_randint(state, 4) ? 6 : 25);
        k = n_randint(state, n_euler_phi(q));
        deflate = n_randint(state, 2);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        /* deflation is not supported for imprimitive principal characters */
        if (q > 1 && dirichlet_conductor_char(G, chi) == 1)
            deflate = 0;

        {
            fmpq_init(x);
            acb_init(s);
            mag_init(mx);
            mag_init(r);
            t = _acb_vec_init(len);
            u = _acb_vec_init(len);

            switch (n_randint(state, 4))
            {
                case 0:
                    fmpq_set_si(x, -8 + (slong) n_randint(state, 12), 1);
                    break;
                case 1:
                    fmpq_set_si(x, -15 + (slong) n_randint(state, 30), 2);
                    break;
                default:
                    fmpq_randtest(x, state, 2 + n_randint(state, 8));
            }

            /* pole */
            if (q == 1 && fmpq_is_one(x))
                deflate = 1;

            /* an inexact s, so that acb_dirichlet_l_jet does not use the
               approximate functional equation */
            acb_set_fmpq(s, x, prec + 20);
            mag_set_ui_2exp_si(arb_radref(acb_realref(s)), 1, -prec - 100);

            acb_dirichlet_l_jet(t, s, G, chi, deflate, len, prec + 20);

            /* the terms are computed in parallel when threads are available */
            flint_set_num_threads(1 + n_randint(state, 4));
            acb_dirichlet_l_fmpq_afe_jet(u, x, (q == 1 && n_randint(state, 2)) ? NULL : G, chi, deflate, len, prec);

            /* the result must not depend on the number of threads */
            if (n_randint(state, 4) == 0)
            {
                acb_ptr v = _acb_vec_init(len);

                flint_set_num_threads(1 + n_randint(state, 4));
                acb_dirichlet_l_fmpq_afe_jet(v, x, G, chi, deflate, len, prec);

                for (j = 0; j < len; j++)
                {
                    if (!acb_equal(u + j, v + j))
                    {
                        flint_printf("FAIL: thread independence\n\n");
                        flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, q, k, prec, len, deflate, j);
                        flint_printf("u = "); acb_printn(u + j, 100, 0); flint_printf("\n\n");
                        flint_printf("v = "); acb_printn(v + j, 100, 0); flint_printf("\n\n");
                        flint_abort();
                    }
                }

                _acb_vec_clear(v, len);
            }

            flint_set_num_threads(1);

            for (j = 0; j < len; j++)
            {
                if (!acb_overlaps(t + j, u + j))
                {
                    flint_printf("FAIL: overlap\n\n");
                    flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, q, k, prec, len, deflate, j);
                    flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                    flint_printf("t = "); acb_printn(t + j, 100, 0); flint_printf("\n\n");
                    flint_printf("u = "); acb_printn(u + j, 100, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }

            /* accuracy relative to the largest coefficient */
            if (prec >= 64)
            {
                for (j = 0; j < len; j++)
                {
                    acb_get_mag(r, t + j);
                    mag_max(mx, mx, r);
                }

                if (!mag_is_zero(mx))
                {
                    mag_mul_2exp_si(mx, mx, -prec + 40);

                    for (j = 0; j < len; j++)
                    {
                        if (mag_cmp(arb_radref(acb_realref(u + j)), mx) > 0 ||
                            mag_cmp(arb_radref(acb_imagref(u + j)), mx) > 0)
                        {
                            flint_printf("FAIL: accuracy\n\n");
                            flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, q, k, prec, len, deflate, j);
                            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                            flint_printf("t = "); acb_printn(t + j, 100, 0); flint_printf("\n\n");
                            flint_printf("u = "); acb_printn(u + j, 100, 0); flint_printf("\n\n");
                            flint_abort();
                        }
                    }
                }
            }

            fmpq_clear(x);
            acb_clear(s);
            mag_clear(mx);
            mag_clear(r);
            _acb_vec_clear(t, len);
            _acb_vec_clear(u, len);
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
    }

    TEST_FUNCTION_END(state);
}
