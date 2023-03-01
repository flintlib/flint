/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"
#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mul_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, z;
        mag_t xb, yb;
        fmpz_t e;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);

        fmpz_init(e);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 100);
        fmpz_randtest(e, state, 100);
        mag_get_fmpr(x, xb);

        mag_mul_2exp_fmpz(yb, xb, e);

        fmpr_mul_2exp_fmpz(y, x, e);

        mag_get_fmpr(z, yb);

        MAG_CHECK_BITS(yb)

        if (!fmpr_equal(z, y))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);

        mag_clear(xb);
        mag_clear(yb);

        fmpz_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
