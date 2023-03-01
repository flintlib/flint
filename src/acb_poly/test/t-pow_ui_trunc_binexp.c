/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("pow_ui_trunc_binexp....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpz_poly */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong zbits1, rbits1, rbits2, trunc;
        ulong e;
        fmpz_poly_t A, B;
        acb_poly_t a, b;

        zbits1 = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        e = n_randint(state, 50);
        trunc = n_randint(state, 40);

        fmpz_poly_init(A);
        fmpz_poly_init(B);

        acb_poly_init(a);
        acb_poly_init(b);

        fmpz_poly_randtest(A, state, 1 + n_randint(state, 10), zbits1);
        fmpz_poly_pow_trunc(B, A, e, trunc);

        acb_poly_set_fmpz_poly(a, A, rbits1);
        acb_poly_pow_ui_trunc_binexp(b, a, e, trunc, rbits2);

        if (!acb_poly_contains_fmpz_poly(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);
            flint_printf("e = %wu\n", e);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("A = "); fmpz_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpz_poly_print(B); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_pow_ui_trunc_binexp(a, a, e, trunc, rbits2);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
