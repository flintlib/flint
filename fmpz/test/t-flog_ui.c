/*
    Copyright (C) 2011 Sebastian Pancratz

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("flog_ui....");
    fflush(stdout);

    

    /* Check correctness */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, x, y;
        slong k;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(x);
        fmpz_init(y);

        while (fmpz_cmp_ui(a, 1) < 0)
            fmpz_randtest(a, state, 200);
        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_set_ui(b, n_randtest(state));

        k = fmpz_flog_ui(a, fmpz_get_ui(b));

        fmpz_pow_ui(x, b, k);
        fmpz_pow_ui(y, b, k + 1);

        result = (fmpz_cmp(x, a) <= 0 && fmpz_cmp(a, y) < 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            flint_printf("k = %wd\n", k);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /* Check correctness:  exact powers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        slong k, l;

        fmpz_init(a);
        fmpz_init(b);

        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_set_ui(b, n_randtest(state));
        l = n_randint(state, 20);
        fmpz_pow_ui(a, b, l);

        k = fmpz_flog_ui(a, fmpz_get_ui(b));

        result = (k == l);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("k = %wd\n", k);
            flint_printf("l = %wd\n", l);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

