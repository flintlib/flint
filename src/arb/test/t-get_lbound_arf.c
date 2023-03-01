/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_lbound_arf....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x;
        arf_t b, b2, b3;
        fmpq_t q;

        arb_init(x);
        arf_init(b);
        arf_init(b2);
        arf_init(b3);
        fmpq_init(q);

        arb_randtest(x, state, 1 + n_randint(state, 200), 10);

        arb_get_lbound_arf(b, x, 2 + n_randint(state, 200));
        arb_get_rand_fmpq(q, state, x, 1 + n_randint(state, 200));
        arf_mul_fmpz(b2, b, fmpq_denref(q), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_fmpz(b3, fmpq_numref(q));

        if (arf_cmp(b2, b3) > 0)
        {
            flint_printf("FAIL (lbound):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("b2 = "); arf_print(b2); flint_printf("\n\n");
            flint_printf("b3 = "); arf_print(b3); flint_printf("\n\n");
            flint_abort();
        }

        arb_get_ubound_arf(b, x, 2 + n_randint(state, 200));
        arb_get_rand_fmpq(q, state, x, 1 + n_randint(state, 200));
        arf_mul_fmpz(b2, b, fmpq_denref(q), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_fmpz(b3, fmpq_numref(q));

        if (arf_cmp(b2, b3) < 0)
        {
            flint_printf("FAIL (ubound):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("b2 = "); arf_print(b2); flint_printf("\n\n");
            flint_printf("b3 = "); arf_print(b3); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        arb_get_lbound_arf(b, x, 2 + n_randint(state, 200));
        arb_get_ubound_arf(b2, x, 2 + n_randint(state, 200));

        if (arf_cmp(b, b2) > 0)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("b2 = "); arf_print(b2); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        fmpq_clear(q);
        arf_clear(b);
        arf_clear(b2);
        arf_clear(b3);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

