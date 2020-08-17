/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2016, 2020 William Hart

    This file is part of FLINT.
    
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "long_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("cmp_si....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        slong y;
        mpq_t X, Y;
        int c1, c2;

        fmpq_init(x);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        y = z_randtest(state);

        fmpq_get_mpq(X, x);
        flint_mpq_set_si(Y, y, 1);

        c1 = fmpq_cmp_si(x, y);
        c2 = mpq_cmp(X, Y);

        if (c1 < 0) c1 = -1;
        if (c1 > 0) c1 = 1;

        if (c2 < 0) c2 = -1;
        if (c2 > 0) c2 = 1;

        if (c1 != c2)
        {
            flint_printf("FAIL\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ny = %wd", y);
            flint_printf("\ncmp(x,y) = %d, cmp(X,Y) = %d\n", c1, c2);
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

