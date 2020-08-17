/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("inv....");
    fflush(stdout);

    /* x = y * z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;
        mpq_t X, Y, Z, YY, ZZ;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        mpq_init(X);
        mpq_init(Y);
        mpq_init(Z);
        mpq_init(YY);
        mpq_init(ZZ);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);
        do { fmpq_randtest(z, state, 200); } while (fmpq_is_zero(z));

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        fmpq_get_mpq(Z, z);

        fmpq_inv(y, z);
        fmpq_inv(x, y);

        if (!fmpq_equal(x, z))
        {
            flint_printf("FAIL: applying inv twice did not give back the input!\n");
            abort();
        }

        if (!fmpq_is_canonical(y) || !fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_inv(Y, Z);
        mpq_inv(X, Z);

        fmpq_get_mpq(YY, y);
        fmpq_get_mpq(ZZ, z);

        if (!mpq_equal(Y, YY) || !mpq_equal(Z, ZZ))
        {
            flint_printf("FAIL: fmpq_inv != mpq_inv\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nz = ");
            fmpq_print(z);
            flint_printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);

        mpq_clear(X);
        mpq_clear(Y);
        mpq_clear(Z);
        mpq_clear(YY);
        mpq_clear(ZZ);
    }

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        mpq_t X, Y;

        fmpq_init(x);
        mpq_init(X);
        mpq_init(Y);

        do { fmpq_randtest(x, state, 200); } while (fmpq_is_zero(x));

        fmpq_get_mpq(X, x);

        fmpq_inv(x, x);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_inv(X, X);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            flint_printf("FAIL: fmpq_mul(x,x,y) != mpq_mul(X,X,Y)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\n");
            abort();
        }

        fmpq_clear(x);

        mpq_clear(X);
        mpq_clear(Y);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
