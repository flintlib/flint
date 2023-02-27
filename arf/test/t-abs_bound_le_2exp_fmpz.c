/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("abs_bound_le_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arf_t x, y;
        fmpz_t b;
        int cmp1, cmp2;

        arf_init(x);
        arf_init(y);
        fmpz_init(b);

        arf_randtest_not_zero(x, state, 2 + n_randint(state, 1000), 100);
        arf_abs_bound_le_2exp_fmpz(b, x);

        arf_one(y);
        arf_mul_2exp_fmpz(y, y, b);

        cmp1 = (arf_cmpabs(x, y) <= 0);

        arf_mul_2exp_si(y, y, -1);

        cmp2 = (arf_cmpabs(y, x) < 0);

        arf_mul_2exp_si(y, y, 1);

        if (!cmp1 || !cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

