/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart

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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    ulong cflags = UWORD(0);

    FLINT_TEST_INIT(state);

    flint_printf("get/set_coeff_fmpq....");
    fflush(stdout);  

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        fmpq_t x, y;
        slong coeff, len;

        fmpq_poly_init(a);
        fmpq_init(x);
        fmpq_init(y);
        len = (slong) (n_randint(state, 100) + 1);

        for (j = 0; j < 50; j++)
        {
            fmpq_randtest(x, state, 200);
            coeff = (slong) n_randint(state, len);
            fmpq_poly_set_coeff_fmpq(a, coeff, x);
            fmpq_poly_get_coeff_fmpq(y, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (fmpq_equal(x, y) && !cflags);
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a     = "), fmpq_poly_debug(a), flint_printf("\n\n");
                flint_printf("coeff = %wd\n\n", coeff);
                flint_printf("len   = %wd\n\n", len);
                flint_printf("cflags = %wu\n\n", cflags);
                flint_printf("x = "), fmpq_print(x), flint_printf("\n");
                flint_printf("y = "), fmpq_print(y), flint_printf("\n");
                abort();
            }
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_poly_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
