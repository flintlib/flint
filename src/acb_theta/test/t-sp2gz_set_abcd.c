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

    flint_printf("sp2gz_set_abcd....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        fmpz_mat_t a, b, c, d;
        fmpz_mat_t m, n;
        slong bits = n_randint(state, 10);

        fmpz_mat_init(a, g, g);
        fmpz_mat_init(b, g, g);
        fmpz_mat_init(c, g, g);
        fmpz_mat_init(d, g, g);
        fmpz_mat_init(m, 2 * g, 2 * g);
        fmpz_mat_init(n, 2 * g, 2 * g);

        sp2gz_randtest(m, state, bits);
        sp2gz_get_a(a, m);
        sp2gz_get_b(b, m);
        sp2gz_get_c(c, m);
        sp2gz_get_d(d, m);
        sp2gz_set_abcd(n, a, b, c, d);

        if (!fmpz_mat_equal(m, n))
        {
            flint_printf("FAIL\n\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\n\n");
            fmpz_mat_print_pretty(n);
            flint_abort();
            flint_printf("\n\n");
        }

        fmpz_mat_clear(a);
        fmpz_mat_clear(b);
        fmpz_mat_clear(c);
        fmpz_mat_clear(d);
        fmpz_mat_clear(m);
        fmpz_mat_clear(n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
