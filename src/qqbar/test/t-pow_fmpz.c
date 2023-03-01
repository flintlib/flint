/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("pow_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check x^m x^n = x^(m+n) */
    for (iter = 0; iter < 100 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, xm, xn, xmxn, xmn;
        fmpz_t m, n, mn;

        qqbar_init(x);
        qqbar_init(xm);
        qqbar_init(xn);
        qqbar_init(xmxn);
        qqbar_init(xmn);
        fmpz_init(m);
        fmpz_init(n);
        fmpz_init(mn);

        if (n_randint(state, 2))
        {
            qqbar_root_of_unity(x, n_randint(state, 100), 1 + n_randint(state, 10));

            fmpz_set_si(m, n_randtest(state));
            fmpz_set_si(n, n_randtest(state));
        }
        else
        {
            do {
                qqbar_randtest(x, state, 4, 10);
                fmpz_randtest(m, state, 5);
                fmpz_randtest(n, state, 5);
            } while (qqbar_is_zero(x) && (fmpz_sgn(m) < 0 || fmpz_sgn(n) < 0));
        }

        fmpz_add(mn, m, n);

        qqbar_pow_fmpz(xm, x, m);
        qqbar_pow_fmpz(xn, x, n);
        qqbar_pow_fmpz(xmn, x, mn);

        qqbar_mul(xmxn, xm, xn);

        if (!qqbar_equal(xmxn, xmn))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("m = "); fmpz_print(m); flint_printf("\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("xm = "); qqbar_print(xm); flint_printf("\n\n");
            flint_printf("xn = "); qqbar_print(xn); flint_printf("\n\n");
            flint_printf("xmxn = "); qqbar_print(xmxn); flint_printf("\n\n");
            flint_printf("xmn = "); qqbar_print(xmn); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(xm);
        qqbar_clear(xn);
        qqbar_clear(xmxn);
        qqbar_clear(xmn);
        fmpz_clear(m);
        fmpz_clear(n);
        fmpz_clear(mn);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

