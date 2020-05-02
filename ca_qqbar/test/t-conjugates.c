/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("conjugates....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, z;
        ca_qqbar_ptr r;
        fmpq_t s;
        slong i, d;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);
        fmpq_init(s);

        ca_qqbar_randtest(x, state, 4, 10);
        d = ca_qqbar_degree(x);

        r = ca_qqbar_vec_init(d);

        ca_qqbar_conjugates(r, x);

        for (i = 0; i < d; i++)
            ca_qqbar_add(y, y, r + i);

        fmpq_set_fmpz_frac(s, CA_QQBAR_COEFFS(x) + d - 1, CA_QQBAR_COEFFS(x) + d);
        fmpq_neg(s, s);

        ca_qqbar_set_fmpq(z, s);

        if (!ca_qqbar_equal(y, z))
        {
            flint_printf("FAIL! %d\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); ca_qqbar_print(z); flint_printf("\n\n");
            for (i = 0; i < d; i++)
            {
                flint_printf("r%wd = ", i); ca_qqbar_print(r + i); flint_printf("\n\n");
            }
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(z);
        ca_qqbar_vec_clear(r, d);
        fmpq_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

