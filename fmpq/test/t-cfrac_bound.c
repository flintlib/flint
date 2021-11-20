/*
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
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("cfrac_bound....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, r;
        fmpz * c;
        slong n, bound;

        fmpq_init(x);
        fmpq_init(r);

        /* Test worst case (quotient of Fibonacci numbers) */
        if (n_randint(state, 50) == 1)
        {
            slong v = 1 + n_randint(state, 1000);
            fmpz_fib_ui(fmpq_numref(x), v + 1);
            fmpz_fib_ui(fmpq_denref(x), v);
        }
        else
        {
            fmpq_randtest(x, state, 1 + n_randint(state, 1000));
        }

        bound = fmpq_cfrac_bound(x);
        c = _fmpz_vec_init(bound + 10);
        n = fmpq_get_cfrac(c, r, x, bound);

        if (n > bound)
        {
            flint_printf("FAIL: length=%wd > bound=%wd\n", n, bound);
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(c, bound + 10);
        fmpq_clear(x);
        fmpq_clear(r);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
