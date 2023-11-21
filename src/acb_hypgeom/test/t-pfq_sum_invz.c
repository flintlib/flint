/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_pfq_sum_invz, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr a, b;
        acb_t z, zinv, s1, s2, t1, t2;
        slong i, p, q, n, prec1, prec2;

        p = n_randint(state, 5);
        q = n_randint(state, 5);
        n = n_randint(state, 300);
        prec1 = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        acb_init(z);
        acb_init(zinv);
        acb_init(s1);
        acb_init(s2);
        acb_init(t1);
        acb_init(t2);

        acb_randtest_special(z, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(s1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(t1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(s2, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(t2, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));

        acb_inv(zinv, z, 1 + n_randint(state, 500));

        a = _acb_vec_init(p);
        b = _acb_vec_init(q);

        for (i = 0; i < p; i++)
            acb_randtest(a + i, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        for (i = 0; i < q; i++)
            acb_randtest(b + i, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));

        if (n_randint(state, 2))
            acb_hypgeom_pfq_sum_invz(s1, t1, a, p, b, q, zinv, z, n, prec1);
        else
            acb_hypgeom_pfq_sum_bs_invz(s1, t1, a, p, b, q, z, n, prec1);

        acb_hypgeom_pfq_sum(s2, t2, a, p, b, q, zinv, n, prec2);

        if (!acb_overlaps(s1, s2) || !acb_overlaps(t1, t2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
            flint_printf("s1 = "); acb_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printd(s2, 15); flint_printf("\n\n");
            flint_printf("t1 = "); acb_printd(t1, 15); flint_printf("\n\n");
            flint_printf("t2 = "); acb_printd(t2, 15); flint_printf("\n\n");
            flint_abort();
        }

        _acb_vec_clear(a, p);
        _acb_vec_clear(b, q);

        acb_clear(z);
        acb_clear(zinv);
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(t1);
        acb_clear(t2);
    }

    TEST_FUNCTION_END(state);
}
