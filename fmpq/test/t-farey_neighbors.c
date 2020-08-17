/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq.h"

int
main(void)
{
    slong i, q, steps;
    FLINT_TEST_INIT(state);

    flint_printf("farey_neighbors....");
    fflush(stdout);

    /* walk from -2 to 2 with a known number of steps */

    for (q = 1; q <= 4 * flint_test_multiplier(); q++)
    {
        fmpq_t s, t, cur, left, right;
        fmpz_t Q;

        fmpq_init(s);
        fmpq_init(t);
        fmpq_init(cur);
        fmpq_init(left);
        fmpq_init(right);

        steps = 0;
        for (i = 1; i <= q; i++)
            steps += 4*n_euler_phi(i);

        fmpz_init_set_ui(Q, q);

        fmpq_set_si(cur, -2, 1);

        for (i = 0; i < steps; i++)
        {
            fmpq_farey_neighbors(left, right, cur, Q);

            fmpq_sub(t, cur, left);
            fmpz_one(fmpq_numref(s));
            fmpz_mul(fmpq_denref(s), fmpq_denref(left), fmpq_denref(cur));
            if (!fmpq_equal(s, t))
            {
                flint_printf("FAIL:\n");
                flint_printf("check left neighbor i = %wd, q = %wd\n", i, q);
                flint_abort();
            }

            fmpq_sub(t, right, cur);
            fmpz_one(fmpq_numref(s));
            fmpz_mul(fmpq_denref(s), fmpq_denref(right), fmpq_denref(cur));
            if (!fmpq_equal(s, t))
            {
                flint_printf("FAIL:\n");
                flint_printf("check right neighbor i = %wd, q = %wd\n", i, q);
                flint_abort();
            }

            fmpq_swap(cur, right);
        }

        fmpq_set_si(right, 2, 1);
        if (!fmpq_equal(cur, right))
        {
            flint_printf("FAIL:\n");
            flint_printf("check end q = %wd\n", q);
            flint_abort();            
        }

        fmpz_clear(Q);

        fmpq_clear(s);
        fmpq_clear(t);
        fmpq_clear(cur);
        fmpq_clear(left);
        fmpq_clear(right);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
