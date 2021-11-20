/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    ulong cflags = UWORD(0);

    FLINT_TEST_INIT(state);

    flint_printf("set_array_mpq....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong j, n = 100;
        fmpq_poly_t f, g;
        mpq_t * a;

        a = (mpq_t *) flint_malloc(n * sizeof(mpq_t));
        for (j = 0; j < n; j++)
            mpq_init(a[j]);

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, n), 200);
        for (j = 0; j < f->length; j++)
            fmpq_poly_get_coeff_mpq(a[j], f, j);

        fmpq_poly_set_array_mpq(g, (const mpq_t *) a, n);

        cflags |= fmpq_poly_is_canonical(g) ? 0 : 1;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpq_poly_debug(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        for (j = 0; j < n; j++)
            mpq_clear(a[j]);
        flint_free(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
