/*
    Copyright (C) 2011 Fredrik Johansson

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
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("set_ui_smod....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, mz;
        mp_limb_t m, r;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(mz);

        do { m = n_randtest(state); } while (m < 2);

        fmpz_set_ui(mz, m);
        fmpz_randtest_mod_signed(a, state, mz);

        r = fmpz_fdiv_ui(a, m);

        fmpz_set_ui_smod(b, r, m);

        if (!fmpz_equal(a, b))
        {
            flint_printf("FAIL:\n");
            flint_printf("a: "); fmpz_print(a); flint_printf("\n");
            flint_printf("m: %wu\n", m);
            flint_printf("r: %wu\n", m);
            flint_printf("b: "); fmpz_print(b); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(mz);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
