/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_CRT, state)
{
    slong i, j;
    int sign;

    fmpz_t input;
    fmpz_t result;
    fmpz_t r1;
    fmpz_t m1;
    fmpz_t mprod;
    fmpz_t r2, m2;

    fmpz_init(input);
    fmpz_init(result);
    fmpz_init(r1);
    fmpz_init(m1);
    fmpz_init(r2);
    fmpz_init(m2);
    fmpz_init(mprod);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong nprimes;

        fmpz_set_ui(m2, n_randtest_prime(state, 0));
        nprimes = 1 + n_randint(state, 4);

        fmpz_set_ui(m1, UWORD(1));
        for (j = 0; j < nprimes; )
        {
            ulong t = n_randtest_prime(state, 0);
            if (t != fmpz_get_ui(m2))
            {
                fmpz_mul_ui(m1, m1, t);
                j++;
            }
        }

        fmpz_mul(mprod, m1, m2);

        sign = n_randint(state, 2);

        if (sign)
            fmpz_randtest_mod_signed(input, state, mprod);
        else
            fmpz_randtest_mod(input, state, mprod);

        fmpz_mod(r1, input, m1);
        fmpz_mod(r2, input, m2);

        if (sign && n_randint(state, 2))
        {
            /* If sign is set, fmpz_CRT allows -m_{1} <= r < m_{1}. Set r < 0.*/
            fmpz_sub(r1, r1, m1);
        }

        fmpz_CRT(result, r1, m1, r2, m2, sign);

        if (!fmpz_equal(result, input) || !_fmpz_is_canonical(result))
        {
            flint_printf("FAIL:\n");
            flint_printf("m1: "); fmpz_print(m1); flint_printf("\n");
            flint_printf("m2: "); fmpz_print(m2); flint_printf("\n");
            flint_printf("m1*m2: "); fmpz_print(mprod); flint_printf("\n");
            flint_printf("input: "); fmpz_print(input); flint_printf("\n");
            flint_printf("r1: "); fmpz_print(r1); flint_printf("\n");
            flint_printf("r2: "); fmpz_print(r2); flint_printf("\n");
            flint_printf("result: "); fmpz_print(result); flint_printf("\n");
            flint_printf("%wd Equalness: %d\n", i, fmpz_equal(result, input));
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(input);
    fmpz_clear(result);
    fmpz_clear(r1);
    fmpz_clear(m1);
    fmpz_clear(r2);
    fmpz_clear(m2);
    fmpz_clear(mprod);

    TEST_FUNCTION_END(state);
}
