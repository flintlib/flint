/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t x;
        acb_t s, t, u;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, k;
        slong prec;

        fmpq_init(x);
        acb_init(s);
        acb_init(t);
        acb_init(u);

        q = 1 + n_randint(state, 50);
        prec = 2 + n_randint(state, 300);
        k = n_randint(state, n_euler_phi(q));

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        fmpq_randtest(x, state, 2 + n_randint(state, 100));
        if (fmpq_get_d(x) < -100.0)
            fmpq_neg(x, x);
        acb_set_fmpq(s, x, prec);

        if (n_randint(state, 2))
            acb_dirichlet_l_hurwitz(t, s, NULL, G, chi, prec);
        else
            acb_dirichlet_l(t, s, G, chi, prec);

        acb_dirichlet_l_fmpq(u, x, G, chi, prec);

        if (!acb_overlaps(t, u))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd\n\n", iter, q, k, prec);
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("s = "); acb_printn(s, 100, 0); flint_printf("\n\n");
            flint_printf("t = "); acb_printn(t, 100, 0); flint_printf("\n\n");
            flint_printf("u = "); acb_printn(u, 100, 0); flint_printf("\n\n");
            acb_sub(t, t, u, prec);
            flint_printf("t - u = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        fmpq_clear(x);
        acb_clear(s);
        acb_clear(t);
        acb_clear(u);
    }

    /* acb_dirichlet_l and acb_dirichlet_l_jet at exact dyadic s, which may
       use the approximate functional equation (for large q also at low
       precision), against the generic algorithms */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t x;
        acb_t s, s2;
        acb_ptr t, u;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, k;
        slong prec, len, n, j;
        int deflate;

        fmpq_init(x);
        acb_init(s);
        acb_init(s2);

        q = 1 + n_randint(state, n_randint(state, 4) ? 30 : 1500);
        prec = 2 + n_randint(state, (q > 300) ? 200 : 600);
        len = 1 + n_randint(state, 5);
        deflate = n_randint(state, 2);
        k = n_randint(state, n_euler_phi(q));

        t = _acb_vec_init(len);
        u = _acb_vec_init(len);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        fmpz_set_si(fmpq_numref(x), (slong) n_randint(state, 40) - 16);
        fmpz_one(fmpq_denref(x));
        fmpz_mul_2exp(fmpq_denref(x), fmpq_denref(x), n_randint(state, 4));
        fmpq_canonicalise(x);
        if (q == 1 && fmpq_is_one(x))
            deflate = 1;

        acb_set_fmpq(s, x, prec);
        /* a ball, for which the generic algorithms are used */
        acb_set(s2, s);
        mag_set_ui_2exp_si(arb_radref(acb_realref(s2)), 1, -prec - 100);

        if (n_randint(state, 2))
        {
            acb_dirichlet_l(u, s, G, chi, prec);
            _acb_dirichlet_l_generic(t, s, G, chi, prec + 20);
            n = 1;
        }
        else
        {
            acb_dirichlet_l_jet(u, s, (q == 1 && n_randint(state, 2)) ? NULL : G, chi, deflate, len, prec);
            acb_dirichlet_l_jet(t, s2, G, chi, deflate, len, prec + 20);
            n = len;
        }

        for (j = 0; j < n; j++)
        {
            if (!acb_overlaps(t + j, u + j))
            {
                flint_printf("FAIL: overlap (dyadic)\n\n");
                flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, q, k, prec, n, deflate, j);
                flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                flint_printf("t = "); acb_printn(t + j, 100, 0); flint_printf("\n\n");
                flint_printf("u = "); acb_printn(u + j, 100, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        fmpq_clear(x);
        acb_clear(s);
        acb_clear(s2);
        _acb_vec_clear(t, len);
        _acb_vec_clear(u, len);
    }

    TEST_FUNCTION_END(state);
}
