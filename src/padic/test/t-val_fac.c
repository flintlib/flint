/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic.h"

TEST_FUNCTION_START(padic_val_fac, state)
{
    int i, result;

    /* Check aliasing for padic_val_fac() */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, p;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_randtest_unsigned(a, state, (flint_bitcnt_t) (1.5 * FLINT_BITS));
        fmpz_set(b, a);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        padic_val_fac(c, b, p);
        padic_val_fac(b, b, p);

        result = fmpz_equal(b, c);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    /* Check correctness for padic_val_fac_ui(), p == 2 */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, p;

        ulong s, t, N;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);

        N = n_randint(state, WORD(1) < 13);
        fmpz_set_ui(p, 2);
        fmpz_fac_ui(a, N);

        s = padic_val_fac_ui_2(N);
        t = fmpz_remove(b, a, p);

        result = (s == t);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("N = %wu\n", N);
            flint_printf("s = %wu\n", s);
            flint_printf("t = %wu\n", t);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Check correctness for padic_val_fac_ui(), any p */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, p;

        ulong s, t, N;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);

        N = n_randint(state, WORD(1) < 13);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, FLINT_BITS - 4), 0));
        fmpz_fac_ui(a, N);

        s = padic_val_fac_ui(N, p);
        t = fmpz_remove(b, a, p);

        result = (s == t);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("N = %wu\n", N);
            flint_printf("s = %wu\n", s);
            flint_printf("t = %wu\n", t);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Compare padic_val_fac_ui() with padic_val_fac() */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, p, t, z;

        ulong s, n;

        fmpz_init(a);
        fmpz_init(p);
        fmpz_init(t);
        fmpz_init(z);

        n = n_randint(state, WORD(1) < 13);
        fmpz_set_ui(z, n);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, FLINT_BITS - 4), 0));
        fmpz_fac_ui(a, n);

        s = padic_val_fac_ui(n, p);
        padic_val_fac(t, z, p);

        result = (fmpz_equal_ui(t, s));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("n = %wu\n", n);
            flint_printf("s = %wu\n", s);
            flint_printf("t = "), fmpz_print(t), flint_printf("\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(p);
        fmpz_clear(t);
        fmpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
