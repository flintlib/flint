/*
     Copyright (C) 2020 Vincent Delecroix
     Copyright (C) 2021 Fredrik Johansson

     This file is part of FLINT.

     FLINT is free software: you can redistribute it and/or modify it under
     the terms of the GNU Lesser General Public License (LGPL) as published
     by the Free Software Foundation; either version 2.1 of the License, or
     (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_add_sub_fmpq, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, s, t;
        fmpq_t c;
        int op, alias;

        fmpq_poly_init(a);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_init(c);

        fmpq_poly_randtest(a, state, 1 + n_randint(state, 4), 1 + n_randint(state, 200));
        fmpq_poly_randtest(s, state, 4, 200);
        fmpq_poly_randtest(t, state, 4, 200);
        fmpq_randtest(c, state, 1 + n_randint(state, 200));

        op = n_randint(state, 3);
        alias = n_randint(state, 3);

        fmpq_poly_set_fmpq(s, c);
        if (op == 0)
            fmpq_poly_add(s, a, s);
        else if (op == 1)
            fmpq_poly_sub(s, a, s);
        else
            fmpq_poly_sub(s, s, a);

        if (alias)
        {
            if (op == 0)
                fmpq_poly_add_fmpq(t, a, c);
            else if (op == 1)
                fmpq_poly_sub_fmpq(t, a, c);
            else
                fmpq_poly_fmpq_sub(t, c, a);
        }
        else
        {
            fmpq_poly_set(t, a);

            if (op == 0)
                fmpq_poly_add_fmpq(t, t, c);
            else if (op == 1)
                fmpq_poly_sub_fmpq(t, t, c);
            else
                fmpq_poly_fmpq_sub(t, c, t);
        }

        if (!fmpq_poly_equal(s, t))
        {
           printf("FAIL (op = %d, alias = %d):\n", op, alias);
           printf("a = "); fmpq_poly_print(a); printf("\n");
           printf("s = "); fmpq_poly_print(s); printf("\n");
           printf("t = "); fmpq_poly_print(t); printf("\n");
           printf("c = "); fmpq_print(c); printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
