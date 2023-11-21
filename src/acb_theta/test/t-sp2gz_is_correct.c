/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sp2gz_is_correct, state)
{
    slong iter;

    /* Test: return 1 on various kinds of symplectic matrices; return 0 on
       non-square of even size */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        fmpz_mat_t a, b, m, n;
        slong r = n_randint(state, 10);
        slong c = n_randint(state, 10);
        slong bits = n_randint(state, 100);

        fmpz_mat_init(a, g, g);
        fmpz_mat_init(b, g, g);
        fmpz_mat_init(m, 2 * g, 2 * g);
        fmpz_mat_init(n, r, c);

        if (iter == 0)
        {
            sp2gz_j(m);
        }
        else if (iter <= sp2gz_nb_fundamental(g))
        {
            sp2gz_fundamental(m, iter - 1);
        }
        else if (iter % 2 == 0)
        {
            fmpz_mat_one(a);
            fmpz_mat_randops(a, state, bits);
            sp2gz_block_diag(m, a);
        }
        else
        {
            fmpz_mat_randtest(a, state, bits);
            fmpz_mat_transpose(b, a);
            fmpz_mat_add(a, a, b);
            sp2gz_trig(m, a);
        }

        if (!sp2gz_is_correct(m))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_one(n);
        if ((r != c || r % 2 == 1) && sp2gz_is_correct(n))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(n);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(a);
        fmpz_mat_clear(b);
        fmpz_mat_clear(m);
        fmpz_mat_clear(n);
    }

    TEST_FUNCTION_END(state);
}

