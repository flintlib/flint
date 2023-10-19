/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2018 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpq.h"

void check_invalid(char * s, int b)
{
    fmpq_t r;
    int err;

    fmpq_init(r);
    err = fmpq_set_str(r, s, b);
    if (!err)
    {
        printf("Got no error with s='%s'\n", s);
        printf("r = "); fmpq_print(r); printf("\n");
        fflush(stdout);
        flint_abort();
    }
    fmpq_clear(r);
}

TEST_FUNCTION_START(fmpq_get_set_str, state)
{
    int i;

    check_invalid("x5/3", 6);
    check_invalid("5x/3", 6);
    check_invalid("5/x3", 6);
    check_invalid("5/3x", 6);

    for (i = 0; i < 100000; i++)
    {
        fmpq_t a, a2;
        mpq_t b;
        int ans, base;
        char *str1, *str2;

        fmpq_init(a);
        fmpq_init(a2);
        mpq_init(b);
        fmpq_randtest(a, state, 200);
        base = (int) (n_randint(state, 31) + 2);

        fmpq_get_mpq(b, a);

        str1 = fmpq_get_str(NULL, base, a);
        str2 = mpq_get_str(NULL, base, b);
        if (strcmp(str1, str2))
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Qd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("str1 = %s\n, str2 = %s\n", str1, str2);
            fflush(stdout);
            flint_abort();
        }

        ans = fmpq_set_str(a2, str1, base);

        if (ans || !fmpq_equal(a, a2))
        {
            flint_printf("FAIL:\n");
            flint_printf("str1 = %s\n", str1);
            flint_printf("base = %d\n", base);
            flint_printf("ans = %d\n", ans);
            fflush(stdout);
            flint_abort();
        }

        flint_free(str1);
        flint_free(str2);

        fmpq_clear(a);
        fmpq_clear(a2);
        mpq_clear(b);
    }

    TEST_FUNCTION_END(state);
}
