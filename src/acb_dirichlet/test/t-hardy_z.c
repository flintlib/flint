/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_hardy_z, state)
{
    slong iter;

    /* test self-consistency */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, s2;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_ptr vec1, vec2;
        slong len1, len2;
        slong prec1, prec2;
        ulong q, k;
        slong i;

        len1 = n_randint(state, 6);
        len2 = n_randint(state, 6);
        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        do {
            q = 1 + n_randint(state, 30);
        } while (q % 4 == 2);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);

        do {
            k = n_randint(state, n_euler_phi(q));
            dirichlet_char_index(chi, G, k);
        } while (dirichlet_conductor_char(G, chi) != q);

        acb_init(s);
        acb_init(s2);
        vec1 = _acb_vec_init(len1);
        vec2 = _acb_vec_init(len2);

        acb_randtest(s, state, 2 + n_randint(state, 200), 2);
        acb_randtest(s2, state, 2 + n_randint(state, 200), 2);
        acb_sub(s2, s2, s2, 200);
        acb_add(s2, s, s2, 200);

        acb_dirichlet_hardy_z(vec1, s, G, chi, len1, prec1);
        acb_dirichlet_hardy_z(vec2, s2, G, chi, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!acb_overlaps(vec1 + i, vec2 + i))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd  q = %wu  k = %wu  i = %wd\n\n", iter, q, k, i);
                flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printn(vec1 + i, 50, 0); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printn(vec2 + i, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (arb_contains_zero(acb_imagref(s)))
        {
            for (i = 0; i < len1; i++)
            {
                if (!arb_contains_zero(acb_imagref(vec1 + i)))
                {
                    flint_printf("FAIL: real 1\n\n");
                    flint_printf("iter = %wd  q = %wu  k = %wu  i = %wd\n\n", iter, q, k, i);
                    flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                    flint_printf("r1 = "); acb_printn(vec1 + i, 50, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        if (arb_contains_zero(acb_imagref(s2)))
        {
            for (i = 0; i < len2; i++)
            {
                if (!arb_contains_zero(acb_imagref(vec2 + i)))
                {
                    flint_printf("FAIL: real 1\n\n");
                    flint_printf("iter = %wd  q = %wu  k = %wu  i = %wd\n\n", iter, q, k, i);
                    flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                    flint_printf("r1 = "); acb_printn(vec2 + i, 50, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        acb_clear(s);
        acb_clear(s2);
        _acb_vec_clear(vec1, len1);
        _acb_vec_clear(vec2, len2);
    }

    TEST_FUNCTION_END(state);
}
