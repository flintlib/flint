/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("sub_si_inline....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpz_t a, c, d;
        slong b;

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest(a, state, 1 + n_randint(state, 200));
        fmpz_randtest(c, state, 1 + n_randint(state, 200));
        fmpz_randtest(d, state, 1 + n_randint(state, 200));
        b = n_randtest(state);

        fmpz_sub_si(c, a, b);
        fmpz_sub_si_inline(d, a, b);

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL\n");
            fmpz_print(a); flint_printf("\n\n");
            flint_printf("%wd", b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_sub_si_inline(a, a, b);
        if (!fmpz_equal(c, a))
        {
            flint_printf("FAIL (aliasing)\n");
            fmpz_print(a); flint_printf("\n\n");
            flint_printf("%wd", b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

