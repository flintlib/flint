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

    flint_printf("floor....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, t, u;
        fmpz_t n, n1;
        int ok;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(t);
        ca_qqbar_init(u);
        fmpz_init(n);
        fmpz_init(n1);

        if (n_randint(state, 2))
        {
            ca_qqbar_randtest(x, state, 4, 200);
        }
        else
        {
            fmpq_t c;
            fmpq_init(c);
            fmpz_randtest(n, state, 400);
            fmpq_randtest(c, state, 400);
            fmpq_add_fmpz(c, c, n);
            fmpz_zero(n);
            ca_qqbar_set_fmpq(x, c);
            ca_qqbar_randtest(y, state, 1, 100);
            ca_qqbar_i(t);
            ca_qqbar_mul(y, y, t);
            ca_qqbar_add(x, x, y);
            fmpq_clear(c);
        }

        ca_qqbar_floor(n, x);

        fmpz_add_ui(n1, n, 1);

        ca_qqbar_set_fmpz(t, n);
        ca_qqbar_set_fmpz(u, n1);

        ok = (ca_qqbar_cmp_re(t, x) <= 0 && ca_qqbar_cmp_re(x, u) < 0);

        if (!ok)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("n1 = "); fmpz_print(n1); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(t);
        ca_qqbar_clear(u);
        fmpz_clear(n);
        fmpz_clear(n1);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

