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

    flint_printf("cos_pi_fmpq_algebraic....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t c1, c2;
        ulong p, q, g;
        slong prec;

        prec = 2 + n_randint(state, 5000);
        q = 1 + n_randint(state, 500);
        p = n_randint(state, q / 2 + 1);

        g = n_gcd(q, p);
        q /= g;
        p /= g;

        arb_init(c1);
        arb_init(c2);

        _arb_cos_pi_fmpq_algebraic(c1, p, q, prec);

        arb_const_pi(c2, prec);
        arb_mul_ui(c2, c2, p, prec);
        arb_div_ui(c2, c2, q, prec);
        arb_cos(c2, c2, prec);

        if (!arb_overlaps(c1, c2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("p/q = %wu/%wu", p, q); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printd(c1, 15); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printd(c2, 15); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(c1) < prec - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("p/q = %wu/%wu", p, q); flint_printf("\n\n");
            flint_printf("prec=%wd eff=%wd\n", prec, arb_rel_accuracy_bits(c1));
            flint_printf("c1 = "); arb_printd(c1, 15); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printd(c2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(c1);
        arb_clear(c2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

