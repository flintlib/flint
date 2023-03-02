/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"


int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t a;
        slong j, k;
        slong rows = n_randint(state, 100);
        slong cols = n_randint(state, 100);

        fmpq_mat_init(a, rows, cols);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                fmpq_zero(fmpq_mat_entry(a, j, k));

        fmpq_mat_clear(a);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
