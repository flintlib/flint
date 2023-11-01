/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_get_mpn, state)
{

    fmpz_t a, b, mmin;
    int i, j, k;
    mp_ptr mpna;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(mmin);

    for (i = 0; i < 100; i++)
    {
        fmpz_set_ui(mmin, 1);
        fmpz_mul_2exp(mmin, mmin, i * FLINT_BITS);

        for (j = 0; j < 100; j++)
        {
            fmpz_set_ui(b, 0);

            k = n_randint(state, FLINT_BITS);
            k += 1;                     /* 1 <= k <= FLINT_BITS */
            k += (i * FLINT_BITS);              /* 2^(i*FLINT_BITS) + 1 <= k <= 2^((i + 1)*FLINT_BITS) */

            fmpz_randtest_unsigned(a, state, k);
            fmpz_add(a, a, mmin);

            k = fmpz_get_mpn(&mpna, a);

            while (k != 0)
            {
                fmpz_add_ui(b, b, mpna[k - 1]);

                if (k > 1)
                    fmpz_mul_2exp(b, b, FLINT_BITS);

                k--;
            }

            if (fmpz_cmp(a, b))
            {
                printf("conversion failed.\nn : ");
                fmpz_print(a);
                printf("\nconverted value : ");
                fmpz_print(b);
                printf("\n");
                fflush(stdout);
                flint_abort();
            }

            flint_free(mpna);
        }
    }

    /* regression test */
    {
#if FLINT_BITS == 32
	slong ksave;
	ulong limb0, limb1;

	fmpz_set_ui(a, 524287);
	fmpz_mul_2exp(a, a, FLINT_BITS);
	fmpz_add_ui(a, a, 4294443520);
        fmpz_zero(b);
        k = fmpz_get_mpn(&mpna, a);

	ksave = k;
	limb0 = mpna[0];
	limb1 = mpna[1];

        while (k != 0)
        {
            fmpz_add_ui(b, b, mpna[k - 1]);

            if (k > 1)
                fmpz_mul_2exp(b, b, FLINT_BITS);

            k--;
        }

        if (fmpz_cmp(a, b))
        {
            flint_printf("limb0 = %wu, limb1 = %wu, k = %wd\n", limb0, limb1, k);
	    printf("conversion failed.\nn : ");
	    fmpz_print(a);
            printf("\nconverted value : ");
            fmpz_print(b);
            printf("\n");
            fflush(stdout);
            flint_abort();
        }

	flint_free(mpna);
#endif
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(mmin);

    TEST_FUNCTION_END(state);
}
