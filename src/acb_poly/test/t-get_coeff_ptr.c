/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("get_coeff_ptr....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000; i++)
    {
        acb_poly_t A;
        acb_t a;
        slong n = n_randint(state, 100);

        acb_poly_init(A);
        acb_poly_randtest(A, state, n_randint(state, 100), 100, 10);
        acb_init(a);

        acb_poly_get_coeff_acb(a, A, n);

        result = n < acb_poly_length(A) ? 
                     acb_equal(a, acb_poly_get_coeff_ptr(A, n)) : 
                     acb_poly_get_coeff_ptr(A, n) == NULL;
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A = "), acb_poly_printd(A, 10), flint_printf("\n\n");
            flint_printf("a = "), acb_print(a), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_abort();
        }

        acb_poly_clear(A);
        acb_clear(a);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

