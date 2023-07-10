/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("sp2gz_is_correct....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        fmpz_mat_t a, b, m;
        slong bits = n_randint(state, 100);

        fmpz_mat_init(a, g, g);
        fmpz_mat_init(b, g, g);
        fmpz_mat_init(m, 2 * g, 2 * g);

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
            flint_printf("FAIL\n\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mat_clear(a);
        fmpz_mat_clear(b);
        fmpz_mat_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

