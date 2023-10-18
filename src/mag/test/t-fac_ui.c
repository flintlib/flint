/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mag.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("fac_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        fmpz_t f;
        mag_t xb;
        ulong n;

        arf_init(x);
        arf_init(y);
        fmpz_init(f);
        mag_init(xb);

        mag_randtest_special(xb, state, 80);
        n = n_randtest(state) % 2000;

        mag_fac_ui(xb, n);
        fmpz_fac_ui(f, n);
        arf_set_fmpz(x, f);
        arf_set_mag(y, xb);

        MAG_CHECK_BITS(xb)

        if (!(arf_cmpabs(y, x) >= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(f);
        mag_clear(xb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
