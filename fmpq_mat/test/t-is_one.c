/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("is_one....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpq_mat_init(A, rows, cols);
        fmpq_mat_one(A);

        if (!fmpq_mat_is_one(A))
        {
            flint_printf("FAIL!\n");
            abort();
        }

        if (rows && cols)
        {
            fmpq_mat_randbits(A, state, 100);
            if (fmpq_mat_is_one(A))
            {
                flint_printf("FAIL!\n");
                abort();
            }
        }

        fmpq_mat_clear(A);
    }



    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
