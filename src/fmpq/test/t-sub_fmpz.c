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
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_sub_fmpz, state)
{
    int i;

    /* x = y + z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        fmpz_t z, one;
        fmpq_t X, Y, Z;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(z);
        fmpz_init_set_ui(one, 1);
        fmpq_init(X);
        fmpq_init(Y);
        fmpq_init(Z);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);
        fmpz_randtest(z, state, 200);

        fmpq_set(X, x);
        fmpq_set(Y, y);
        fmpq_set_fmpz_frac(Z, z, one);

        fmpq_sub_fmpz(x, y, z);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_sub(X, Y, Z);

        if (!fmpq_equal(X, x))
        {
            flint_printf("FAIL: fmpq_add(x,y,z) != mpq_add(X,Y,Z)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nz = ");
            fmpz_print(z);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(z);
        fmpz_clear(one);

        fmpq_clear(X);
        fmpq_clear(Y);
        fmpq_clear(Z);
    }

    /* x = x + y */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        fmpz_t y, one;
        fmpq_t X, Y;

        fmpq_init(x);
        fmpz_init(y);
        fmpz_init(one);
        fmpq_init(X);
        fmpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpz_randtest(y, state, 200);
        fmpz_init_set_ui(one, 1);

        fmpq_set(X, x);
        fmpq_set_fmpz_frac(Y, y, one);

        fmpq_sub_fmpz(x, x, y);
        fmpq_sub(X, X, Y);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpq_equal(X, x))
        {
            flint_printf("FAIL: fmpq_add(x,x,y) != mpq_add(X,X,Y)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpz_print(y);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpz_clear(y);
        fmpz_clear(one);

        fmpq_clear(X);
        fmpq_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
