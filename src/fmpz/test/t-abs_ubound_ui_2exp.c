/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_abs_ubound_ui_2exp, state)
{
    slong iter;
    int result;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t x, y;
        slong bits, yexp;
        slong exp;
        mp_limb_t man;

        fmpz_init(x);
        fmpz_init(y);

        fmpz_randtest_not_zero(x, state, 1 + n_randint(state, 400));

        bits = 1 + n_randint(state, FLINT_BITS - 1);

        /* compute an exactly rounded mantissa */
        fmpz_abs(y, x);

        if (fmpz_is_zero(y))
        {
            yexp = 0;
        }
        else
        {
            yexp = fmpz_bits(y) - bits;

            if (yexp >= 0)
            {
                fmpz_cdiv_q_2exp(y, y, yexp);
                if (fmpz_bits(y) == bits + 1)
                {
                    fmpz_tdiv_q_2exp(y, y, 1);
                    yexp--;
                }
            }
            else
            {
                fmpz_mul_2exp(y, y, -yexp);
            }
        }

        man = fmpz_abs_ubound_ui_2exp(&exp, x, bits);

        if (FLINT_BIT_COUNT(man) != bits)
        {
            flint_printf("wrong number of bits!\n");
            flint_printf("bits = %wd, count = %u\n\n", bits, FLINT_BIT_COUNT(man));
            flint_printf("x = "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("bits(x) = %wd\n\n", fmpz_bits(x));
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_printf("yexp = %wd\n\n", yexp);
            flint_printf("man = %wu, exp = %wd\n", man, exp);
            fflush(stdout);
            flint_abort();
        }

        /* ok if equal */
        result = (fmpz_cmp_ui(y, man) == 0);

        /* ok if mantissa is 1 larger */
        if (!result)
        {
            result = ((exp == yexp) && (fmpz_cmp_ui(y, man - 1) == 0));
        }

        /* ok if the exact mantissa is 2^r-1 and overflow to 2^r happened */
        if (!result)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_ui(t, man);
            fmpz_mul_ui(t, t, 2);
            fmpz_sub_ui(t, t, 1);
            result = (exp == yexp + 1) && fmpz_equal(t, y);
            fmpz_clear(t);
        }

        if (!result)
        {
            flint_printf("different from exact ceiling division\n");
            flint_printf("bits = %wd\n\n", bits);
            flint_printf("x = "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("bits(x) = %wd\n\n", fmpz_bits(x));
            flint_printf("y = "); fmpz_print(y); flint_printf(", yexp = %wd\n\n", yexp);
            flint_printf("man = %wu, exp = %wd\n", man, exp);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
