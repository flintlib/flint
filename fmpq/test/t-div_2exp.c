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
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("div_2exp....");
    fflush(stdout);

    /* x = y * 2^exp */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;
        flint_bitcnt_t c;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_numref(y), fmpq_numref(y), n_randint(state, 200));
        else if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), n_randint(state, 200));
        fmpq_canonicalise(y);

        c = n_randint(state, 200);
        fmpq_div_2exp(x, y, c);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical!\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nc = %wu\n", c);
            abort();
        }

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        mpq_div_2exp(X, Y, c);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            flint_printf("FAIL: fmpq_div_2exp(x,y,c) != mpq_div_2exp(X,Y,c)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        mpq_clear(X);
        mpq_clear(Y);
    }

    /* y = y * 2^exp */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        flint_bitcnt_t c;

        fmpq_init(x);
        fmpq_init(y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_numref(y), fmpq_numref(y), n_randint(state, 200));
        else if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), n_randint(state, 200));
        fmpq_canonicalise(y);

        c = n_randint(state, 200);
        fmpq_div_2exp(x, y, c);
        fmpq_div_2exp(y, y, c);

        if (!fmpq_is_canonical(y))
        {
            flint_printf("FAIL: result not canonical!\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nc = %wu\n", c);
            abort();
        }

        if (!fmpq_equal(x, y))
        {
            flint_printf("FAIL: fmpq_div_2exp(x,y,c) != fmpq_div_2exp(y,y,c)\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = ");
            fmpq_print(y);
            flint_printf("\nc = %wu\n", c);
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
