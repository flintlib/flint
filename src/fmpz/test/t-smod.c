/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_smod, state)
{
    int i;

    for (i = 0; i < 20000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, e;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        if (n_randint(state, 10) == 0)
        {
            fmpz_fdiv_q_ui(a, b, 2);
            if (n_randint(state, 2))
                fmpz_add(a, b, b);
            if (n_randint(state, 2))
                fmpz_sub(a, b, b);
            if (n_randint(state, 2))
                fmpz_add(a, b, b);
        }

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpz_smod(c, a, b);
        }
        else if (aliasing == 1)
        {
            fmpz_set(a, b);
            fmpz_smod(c, a, a);
        }
        else if (aliasing == 2)
        {
            fmpz_set(c, a);
            fmpz_smod(c, c, b);
        }
        else
        {
            fmpz_set(c, b);
            fmpz_smod(c, a, c);
        }

        fmpz_sub(d, a, c);
        fmpz_mod(d, d, b);
        if (!fmpz_is_zero(d) || !_fmpz_is_canonical(c))
        {
            flint_printf("FAIL: check b|(smod(a,b) - a)\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_abs(e, b);
        fmpz_mul_2exp(d, c, 1);
        if (fmpz_cmp(d, e) > 0)
        {
            flint_printf("FAIL: check 2*smod(a,b) <= |b|\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_neg(e, e);
        if (fmpz_cmp(e, d) >= 0)
        {
            flint_printf("FAIL: check -|b| < 2*smod(a,b)\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
