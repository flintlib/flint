/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_char_dot, state)
{
    slong iter;

    /* Test: various dots are the same */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong prec = 100;
        ulong a, b;
        slong * n;
        acb_ptr v, w;
        fmpz_t m;
        slong x1, x2;
        acb_t x3, x4;
        int res;

        n = flint_malloc(g * sizeof(slong));
        v = _acb_vec_init(g);
        w = _acb_vec_init(g);
        a = n_randint(state, 1 << g);
        b = n_randint(state, 1 << g);
        acb_init(x3);
        acb_init(x4);
        fmpz_init(m);

        x1 = acb_theta_char_dot(a, b, g);
        acb_theta_char_get_slong(n, b, g);
        x2 = acb_theta_char_dot_slong(a, n, g);
        acb_theta_char_get_acb(v, b, g);
        acb_theta_char_dot_acb(x3, a, v, g, prec);
        acb_theta_char_get_acb(w, a, g);
        acb_dot(x4, NULL, 0, v, 1, w, 1, g, prec);

        if (x1 != x2
            || !acb_overlaps(x3, x4))
        {
            flint_printf("FAIL\n");
            flint_printf("x1 = %wd, x2 = %wd, x3, x4:\n", x1, x2);
            acb_printd(x3, 10);
            flint_printf("\n");
            acb_printd(x4, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mul_2exp_si(x3, x3, 2);
        res = acb_get_unique_fmpz(m, x3);
        fmpz_sub_si(m, m, x1);
        if (!res || !fmpz_divisible_si(m, 4))
        {
            flint_printf("FAIL (mod 4)\n");
            flint_abort();
        }

        flint_free(n);
        _acb_vec_clear(v, g);
        _acb_vec_clear(w, g);
        acb_clear(x3);
        acb_clear(x4);
        fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
