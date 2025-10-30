/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"

typedef struct
{
    int base;
    const char* input;
    int rv;
    slong result;
}
testcases_t;

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

    /* Test specific cases to ensure proper edge-case handling */
    testcases_t cases[] =
    {
        {10, "", -1, 0},
        { 0, "", -1, 0},
        {10, "0", 0, 0},
        { 0, "0", 0, 0},
        {10, "-", -1, 0},
        { 0, "-", -1, 0},
        {10, "-0", 0, 0},
        { 0, "-0", 0, 0},
        {10, "- 2", 0, -2},
        { 0, "- 2", 0, -2},
        {10, " -2", 0, -2},
        { 0, " -2", 0, -2},
        {10, " - 2", 0, -2},
        { 0, " - 2", 0, -2},
        {10, "-  2", 0, -2},
        { 0, "-  2", 0, -2},
        {10, "- 2 ", 0, -2},
        { 0, "- 2 ", 0, -2},
        {10, "--2", -1, 0},
        { 0, "--2", -1, 0},
        {10, "1 2 3", 0, 123},
        { 0, "1 2 3", 0, 123},
        {10, "1-2", -1, 0},
        { 0, "1-2", -1, 0},
        {10, "f", -1, 0},
        { 0, "f", -1, 0},
        { 0, "0x0", 0, 0},
        { 0, "-0x0", 0, 0},
        { 0, "0xf", 0, 15},
        { 0, "-0xf", 0, -15},
        { 0, "--0xf", -1, 0},
        { 0, "0x f", 0, 15},
        { 0, "-0x f", 0, -15},
        { 0, "- 0x f", 0, -15},
        { 0, "0 xf", -1, 0},
        {16, "0x10", -1, 0},
        { 0, "0x10", 0, 16},
        { 2, "0b10", -1, 0},
        { 0, "0b10", 0,  2},
        { 8, "010", 0, 8},
        { 0, "010", 0, 8},
        {10, "010", 0, 10},
    };
    for (i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
    {
        fmpz_t a;
        int ret;

        fmpz_init(a);
        ret = fmpz_set_str(a, cases[i].input, cases[i].base);

        if (ret != cases[i].rv
                || (ret == 0 && fmpz_get_si(a) != cases[i].result))
        {
            flint_printf("FAIL:\n");
            flint_printf("case = %d\n", i);
            flint_printf("base = %d, str = %s\n", cases[i].base, cases[i].input);
            flint_printf("rv = %d, expected = %d\n", ret, cases[i].rv);
            flint_printf("result = "); fmpz_print(a);
            flint_printf(", expected = %ld\n", cases[i].result);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
