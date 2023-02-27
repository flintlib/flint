/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("primorial....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        ulong n, m;
        fmpz_t y;
        arb_t x;
        slong prec;

        arb_init(x);
        fmpz_init(y);

        if (n_randint(state, 10) == 0)
            n = n_randint(state, 5000);
        else
            n = n_randint(state, 500);
        prec = 2 + n_randint(state, 500);

        arb_primorial_ui(x, n, prec);
        fmpz_primorial(y, n);

        if (!arb_contains_fmpz(x, y))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("x = "); arb_printd(x, 100); flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_abort();
        }

        m = n_randint(state, 500);
        if (m == 0)
            n = 0;
        else
            n = n_nth_prime(m);
        prec = 2 + n_randint(state, 500);

        arb_primorial_nth_ui(x, m, prec);
        fmpz_primorial(y, n);

        if (!arb_contains_fmpz(x, y))
        {
            flint_printf("FAIL: containment (2)\n\n");
            flint_printf("m = %wu\n", m);
            flint_printf("x = "); arb_printd(x, 100); flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        fmpz_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

