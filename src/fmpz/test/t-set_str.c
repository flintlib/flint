/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_set_str, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t b;
        char * str;
        int ret1, ret2;
        int base;
        slong i;

        fmpz_init(a);
        fmpz_init(c);
        mpz_init(b);

        base = n_randint(state, 20);
        if (base == 1)
            base = 10;

        fmpz_randtest(a, state, 100);
        str = fmpz_get_str(NULL, base ? base : 10, a);

        /* test inner whitespace */
        i = n_randint(state, strlen(str));
        /* work around GMP bug: "- 123" is not accepted */
        if (n_randint(state, 10) == 0 && i != 0 && str[i - 1] != '-')
            str[i] = ' ';

        /* test error handling */
        if (n_randint(state, 10) == 0)
            str[n_randint(state, strlen(str))] = '?';

        if (n_randint(state, 100) == 0)
            str[0] = '\0';

        fmpz_zero(a);
        ret1 = fmpz_set_str(a, str, base);
        ret2 = mpz_set_str(b, str, base);
        fmpz_set_mpz(c, b);

        if (ret1 != ret2 || (ret1 == 0 && !fmpz_equal(a, c)) || !_fmpz_is_canonical(a))
        {
            flint_printf("FAIL:\n");
            flint_printf("base = %d\n", base);
            flint_printf("str = %s\n", str);
            flint_printf("ret1 = %d, ret2 = %d\n", ret1, ret2);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(str);

        fmpz_clear(a);
        fmpz_clear(c);
        mpz_clear(b);
    }

    /* test binary splitting code */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t b;
        char * str;
        slong i, n;
        int ret1, ret2;

        flint_set_num_threads(1 + n_randint(state, 5));

        fmpz_init(a);
        fmpz_init(c);
        mpz_init(b);

        n = 20000 + n_randint(state, 100000);
        str = flint_malloc(n + 1);

        for (i = 0; i < n; i++)
            str[i] = n_randint(state, 10) + '0';

        str[n] = '\0';

        ret1 = fmpz_set_str(a, str, 10);
        ret2 = mpz_set_str(b, str, 10);
        fmpz_set_mpz(c, b);

        if (ret1 != ret2 || !fmpz_equal(a, c) || !_fmpz_is_canonical(a))
        {
            flint_printf("FAIL:\n");
            flint_printf("str = %s\n", str);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(str);

        fmpz_clear(a);
        fmpz_clear(c);
        mpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
