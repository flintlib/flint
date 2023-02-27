/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "fexpr.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("set_fmpz...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * calcium_test_multiplier(); iter++)
    {
        fmpz_t a, b;
        fexpr_t x;
        slong i;

        fmpz_init(a);
        fmpz_init(b);
        fexpr_init(x);

        for (i = 0; i < 5; i++)
        {
            fmpz_randtest(a, state, 300);
            fexpr_set_fmpz(x, a);
            fexpr_get_fmpz(b, x);

            if (!fmpz_equal(a, b))
            {
                flint_printf("FAIL\n\n");
                flint_printf("a = "); fmpz_print(a); printf("\n");
                flint_printf("b = "); fmpz_print(b); printf("\n");
                flint_printf("x = "); fexpr_print(x); printf("\n");
                flint_abort();
            }
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fexpr_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
