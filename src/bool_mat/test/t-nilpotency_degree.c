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

    flint_printf("nilpotency_degree....");
    fflush(stdout);

    flint_randinit(state);

    /* empty matrix */
    {
        bool_mat_t A;
        bool_mat_init(A, 0, 0);
        if (bool_mat_nilpotency_degree(A) != 0)
        {
            flint_printf("FAIL (empty)\n");
            flint_abort();
        }
        bool_mat_clear(A);
    }

    /* check nilpotency degree by looking at each power of a matrix */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong m, i;
        bool_mat_t A, B;
        slong degree;

        m = n_randint(state, 10) + 1;

        bool_mat_init(A, m, m);
        bool_mat_init(B, m, m);

        bool_mat_randtest(A, state);
        bool_mat_one(B);

        for (degree = -1, i = 0; i < m && degree < 0; i++)
        {
            bool_mat_mul(B, B, A);
            if (!bool_mat_any(B))
                degree = i+1;
        }

        if (bool_mat_nilpotency_degree(A) != degree)
        {
            flint_printf("FAIL (matrix powers)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n");
            flint_printf("B = \n"); bool_mat_print(B); flint_printf("\n");
            flint_printf("i = %wd\n", i);
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
