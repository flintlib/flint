/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("randtest....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x;
        slong deg, bits;
        int real;

        qqbar_init(x);

        real = n_randint(state, 3);
        deg = 1 + n_randint(state, 4) + (real == 2);
        bits = 1 + n_randint(state, 4);

        if (real == 0)
            qqbar_randtest(x, state, deg, bits);
        else if (real == 1)
            qqbar_randtest_real(x, state, deg, bits);
        else if (real == 2)
            qqbar_randtest_nonreal(x, state, deg, bits);

        if (qqbar_degree(x) > deg || qqbar_height_bits(x) > bits ||
            (real == 1 && !qqbar_is_real(x)) || (real == 2 && qqbar_is_real(x)))
        {
            flint_printf("FAIL\n");
            flint_printf("deg = %wd, bits = %wd, real = %d\n", deg, bits, real);
            qqbar_print(x);
            flint_printf("\n");
            flint_abort();
        }

        qqbar_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
