/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("set_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits, res;
        fmpr_t x, z;
        fmpq_t y;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(z);
        fmpq_init(y);

        fmpr_randtest(x, state, bits, 10);
        fmpr_randtest(z, state, bits, 10);

        fmpr_get_fmpq(y, x);
        res = fmpr_set_fmpq(z, y, bits, FMPR_RND_DOWN);

        if (!fmpr_equal(x, z) || !fmpr_check_ulp(z, res, bits))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits: %wd\n", bits);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpq_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
