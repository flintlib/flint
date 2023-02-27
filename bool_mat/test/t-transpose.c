/*
    Copyright (C) 2016 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("transpose....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong m, n;
        bool_mat_t a, b;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bool_mat_init(a, m, n);
        bool_mat_init(b, n, m);

        bool_mat_randtest(a, state);
        bool_mat_randtest(b, state);

        bool_mat_transpose(b, a);

        /* involution */
        {
            bool_mat_t c;
            bool_mat_init(c, m, n);
            bool_mat_randtest(c, state);
            bool_mat_transpose(c, b);
            if (!bool_mat_equal(c, a))
            {
                flint_printf("FAIL (involution)\n");
                flint_printf("m = %wd, n = %wd\n", m, n);
                flint_abort();
            }
            bool_mat_clear(c);
        }

        /* aliasing */
        if (bool_mat_is_square(a))
        {
            bool_mat_transpose(a, a);

            if (!bool_mat_equal(a, b))
            {
                flint_printf("FAIL (aliasing)\n");
                flint_abort();
            }
        }

        bool_mat_clear(a);
        bool_mat_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
