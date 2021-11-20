/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i;
    slong max_limbs = 100;
    ulong * limbs;

    FLINT_TEST_INIT(state);

    flint_printf("get_set_ui_array....");
    fflush(stdout);

    limbs = (ulong *) flint_malloc(max_limbs*sizeof(ulong));

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        slong j;
        slong limb_count1, limb_count2, limb_count3;
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        limb_count1 = n_randint(state, max_limbs) + 1;
        limb_count2 = n_randint(state, 2*max_limbs);
        limb_count3 = n_randint(state, 2*max_limbs);

        fmpz_randtest(a, state, limb_count1*FLINT_BITS);
        fmpz_abs(a, a);
        fmpz_randtest(b, state, limb_count2*FLINT_BITS);
        fmpz_randtest(c, state, limb_count3*FLINT_BITS);

        fmpz_get_ui_array(limbs, limb_count1, a);
        fmpz_set_ui_array(b, limbs, limb_count1);
        fmpz_zero(c);
        for (j = limb_count1 - 1; j >= 0; j--)
        {
            fmpz_mul_2exp(c, c, FLINT_BITS);
            fmpz_add_ui(c, c, limbs[j]);
        }

        if (!fmpz_equal(a, b))
        {
            flint_printf("FAIL:\n");
            flint_printf("Check get and set are inverse\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_equal(a, c))
        {
            flint_printf("FAIL:\n");
            flint_printf("Check limbs are accurate\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    flint_free(limbs);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
