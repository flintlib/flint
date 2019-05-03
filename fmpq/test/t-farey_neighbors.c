/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("farey_neighbors....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t a, s, t, left, right;

        fmpq_init(a);
        fmpq_init(s);
        fmpq_init(t);
        fmpq_init(left);
        fmpq_init(right);

        fmpq_randtest(a, state, 999);

        if (fmpq_farey_neighbors(left, right, a))
        {
            fmpq_sub(t, a, left);
            fmpz_one(fmpq_numref(s));
            fmpz_mul(fmpq_denref(s), fmpq_denref(left), fmpq_denref(a));
            if (fmpz_cmp(fmpq_denref(left), fmpq_denref(a)) >= 0
                || !fmpq_equal(s, t))
            {
                flint_printf("FAIL:\n");
                flint_printf("check left neighbor i = %wd\n", i);
                flint_abort();
            }

            fmpq_sub(t, right, a);
            fmpz_one(fmpq_numref(s));
            fmpz_mul(fmpq_denref(s), fmpq_denref(right), fmpq_denref(a));
            if (fmpz_cmp(fmpq_denref(right), fmpq_denref(a)) >= 0
                || !fmpq_equal(s, t))
            {
                flint_printf("FAIL:\n");
                flint_printf("check right neighbor i = %wd\n", i);
                flint_abort();
            }
        }
        else
        {
            if (!fmpz_is_one(fmpq_denref(a)))
            {
                flint_printf("FAIL:\n");
                flint_printf("check neighbors could be computed i = %wd\n", i);
                flint_abort();
            }
        }

        fmpq_clear(a);
        fmpq_clear(s);
        fmpq_clear(t);
        fmpq_clear(left);
        fmpq_clear(right);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
