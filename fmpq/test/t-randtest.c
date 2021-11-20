/*
    Copyright (C) 2021 Fredrik Johansson

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
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("randtest....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        slong bits;

        fmpq_init(x);
        bits = 1 + n_randint(state, 100);

        fmpq_randtest(x, state, bits);

        if (!fmpq_is_canonical(x) || fmpz_bits(fmpq_numref(x)) > bits || fmpz_bits(fmpq_denref(x)) > bits)
        {
            flint_printf("FAIL\n");
            flint_printf("x: "); fmpq_print(x); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
