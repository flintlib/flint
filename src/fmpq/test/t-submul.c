/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

void mpq_submul(mpq_t x, mpq_t y, mpq_t z)
{
    mpq_t t;
    mpq_init(t);
    mpq_mul(t, y, z);
    mpq_sub(x, x, t);
    mpq_clear(t);
}

TEST_FUNCTION_START(fmpq_submul, state)
{
    int i;

    /* x -= y * z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;
        mpq_t X, Y, Z;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        mpq_init(X);
        mpq_init(Y);
        mpq_init(Z);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);
        fmpq_randtest(z, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        fmpq_get_mpq(Z, z);

        fmpq_submul(x, y, z);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            fflush(stdout);
            flint_abort();
        }

        mpq_submul(X, Y, Z);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            flint_printf("FAIL: fmpq_submul(x,y,z) != mpq_submul(X,Y,Z)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nz = ");
            fmpq_print(z);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);

        mpq_clear(X);
        mpq_clear(Y);
        mpq_clear(Z);
    }

    /* x -= x * y */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);

        fmpq_submul(x, x, y);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            fflush(stdout);
            flint_abort();
        }

        mpq_submul(X, X, Y);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            flint_printf("FAIL: fmpq_submul(x,x,y) != mpq_submul(X,X,Y)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);

        mpq_clear(X);
        mpq_clear(Y);
    }

    /* x -= y * x */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);

        fmpq_submul(x, y, x);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            fflush(stdout);
            flint_abort();
        }

        mpq_submul(X, Y, X);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            flint_printf("FAIL: fmpq_submul(x,y,x) != mpq_submul(X,Y,X)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);

        mpq_clear(X);
        mpq_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
