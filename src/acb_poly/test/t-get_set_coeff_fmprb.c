/*
    Copyright (C) 2009 William Hart
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
    int i, j, result;
    flint_rand_t state;

    flint_printf("get/set_coeff_acb....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        acb_poly_t a;
        acb_t x1, x2;
        slong coeff, len;

        acb_poly_init(a);
        acb_init(x1);
        acb_init(x2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 100; j++)
        {
            acb_randtest(x1, state, 2 + n_randint(state, 200), 10);
            coeff = n_randint(state, len);
            acb_poly_set_coeff_acb(a, coeff, x1);
            acb_poly_get_coeff_acb(x2, a, coeff);

            result = (acb_equal(x1, x2));
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("x1 = "), acb_print(x1), flint_printf("\n");
                flint_printf("x2 = "), acb_print(x2), flint_printf("\n");
                flint_printf("coeff = %wd, length = %wd\n", coeff, len);
                flint_abort();
            }
        }

        acb_clear(x1);
        acb_clear(x2);
        acb_poly_clear(a);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

