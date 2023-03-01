/*
    Copyright (C) 2013 Fredrik Johansson

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

    flint_printf("sin_pi_fmpq_algebraic....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t s1, s2;
        ulong p, q, g;
        slong prec;

        prec = 2 + n_randint(state, 5000);
        q = 1 + n_randint(state, 500);
        p = n_randint(state, q / 2 + 1);

        g = n_gcd(q, p);
        q /= g;
        p /= g;

        arb_init(s1);
        arb_init(s2);

        _arb_sin_pi_fmpq_algebraic(s1, p, q, prec);

        arb_const_pi(s2, prec);
        arb_mul_ui(s2, s2, p, prec);
        arb_div_ui(s2, s2, q, prec);
        arb_sin(s2, s2, prec);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("p/q = %wu/%wu", p, q); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printd(s2, 15); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(s1) < prec - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("p/q = %wu/%wu", p, q); flint_printf("\n\n");
            flint_printf("prec=%wd eff=%wd\n", prec, arb_rel_accuracy_bits(s1));
            flint_printf("s1 = "); arb_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printd(s2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

