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

    flint_printf("exp_pi_i....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x;
        acb_t z, w;
        fmpq_t t;
        slong p;
        ulong q;
        slong prec;

        qqbar_init(x);
        acb_init(z);
        acb_init(w);
        fmpq_init(t);

        q = 1 + n_randint(state, 20);
        p = n_randint(state, 1000);
        p -= 500;
        prec = 2 + n_randint(state, 1000);

        fmpq_set_si(t, p, q);
        arb_sin_cos_pi_fmpq(acb_imagref(w), acb_realref(w), t, prec);

        qqbar_exp_pi_i(x, p, q);
        qqbar_get_acb(z, x, prec);

        if (!acb_overlaps(z, w))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 200, 0); flint_printf("\n\n");
            flint_printf("w = "); acb_printn(w, 200, 0); flint_printf("\n\n");
            flint_printf("p, q = %wd %wu\n\n", p, q);
            flint_abort();
        }

        qqbar_clear(x);
        acb_clear(z);
        acb_clear(w);
        fmpq_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

