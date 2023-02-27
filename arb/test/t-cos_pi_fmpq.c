/*
    Copyright (C) 2013 Fredrik Johansson

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

    flint_printf("cos_pi_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t c1, c2;
        fmpq_t x;
        slong prec;

        prec = 2 + n_randint(state, 5000);

        arb_init(c1);
        arb_init(c2);
        fmpq_init(x);

        fmpq_randtest(x, state, 1 + n_randint(state, 200));

        arb_cos_pi_fmpq(c1, x, prec);

        arb_const_pi(c2, prec);
        arb_mul_fmpz(c2, c2, fmpq_numref(x), prec);
        arb_div_fmpz(c2, c2, fmpq_denref(x), prec);
        arb_cos(c2, c2, prec);

        if (!arb_overlaps(c1, c2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printd(c1, 15); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printd(c2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(c1);
        arb_clear(c2);
        fmpq_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

