/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2018 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

void check_invalid(char * s)
{
    fmpq_poly_t p;
    int err;

    fmpq_poly_init(p);
    err = fmpq_poly_set_str(p, s);
    if (!err)
    {
        printf("Got no error with s='%s'\n", s);
        printf("p = "); fmpq_poly_print(p); printf("\n");
        fflush(stdout);
        flint_abort();
    }
    fmpq_poly_clear(p);
}

TEST_FUNCTION_START(fmpq_poly_get_set_str, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* badly formatted input */
    check_invalid("");
    check_invalid("1");
    check_invalid("x");
    check_invalid("-");
    check_invalid("-1");
    check_invalid("-1  0");
    check_invalid("2 X1 0");
    check_invalid("3 X-2 0 1");
    check_invalid("2 -1 0 1Y");
    check_invalid("2 -1 0 1");
    check_invalid("3   -1 0 1 ");
    check_invalid("3  -1 0 1 ");
    check_invalid("3  -1 0  1");
    check_invalid("3  -1  0 1");

    /* wrong length */
    check_invalid("0  0");
    check_invalid("2  -1 0 1");
    check_invalid("4  0 0");

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int ans;
        char * str;
        fmpq_poly_t f, g;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);

        str = fmpq_poly_get_str(f);
        ans = fmpq_poly_set_str(g, str);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        result = (ans == 0 && fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f      = "), fmpq_poly_debug(f), flint_printf("\n\n");
            flint_printf("g      = "), fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("str    = %s\n\n", str);
            flint_printf("ans    = %d\n\n", ans);
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        flint_free(str);
    }

    TEST_FUNCTION_END(state);
}
