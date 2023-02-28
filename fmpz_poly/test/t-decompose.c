/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_poly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("decompose....");
    fflush(stdout);

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h, G, H, F;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(F);
        fmpz_poly_init(G);
        fmpz_poly_init(H);

        for (j = 0; j < 10; j++)
        {
            fmpz_poly_randtest(g, state, n_randint(state, 35), 80);
            fmpz_poly_randtest(h, state, n_randint(state, 15), 50);

            if (fmpz_poly_degree(g) > 1 && fmpz_poly_degree(h) > 1)
            {
                fmpz_poly_compose(f, g, h);
                if (!fmpz_poly_decompose(G, fmpz_poly_degree(g),
                                         H, fmpz_poly_degree(h), f))
                {
                    flint_printf("FAIL:\n");
                    printf("check decompose success\n");
                    printf("g: "); fmpz_poly_print_pretty(g, "x"); printf("\n");
                    printf("h: "); fmpz_poly_print_pretty(h, "x"); printf("\n");
                    flint_abort();
                }

                fmpz_poly_compose(F, G, H);
                if (!fmpz_poly_equal(F, f))
                {
                    flint_printf("FAIL:\n");
                    printf("check matching decomposition\n");
                    flint_abort();
                }
            }
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(F);
        fmpz_poly_clear(G);
        fmpz_poly_clear(H);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
