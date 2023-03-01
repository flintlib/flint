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

    flint_printf("zeta_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 400 * arb_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        int deflate;
        acb_poly_t S, A, B, C, D, E, F;
        acb_t a, a1;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 25);
        n1 = 1 + n_randint(state, 25);
        n2 = 1 + n_randint(state, 25);

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);
        acb_poly_init(D);
        acb_poly_init(E);
        acb_poly_init(F);
        acb_init(a);
        acb_init(a1);

        deflate = n_randint(state, 2);

        acb_poly_randtest(S, state, m, bits1, 3);
        arb_randtest_precise(acb_realref(a), state, bits1, 3);
        arb_randtest_precise(acb_imagref(a), state, bits1, 3);
        acb_poly_set_coeff_acb(S, 0, a);

        if (n_randint(state, 2))
            acb_randtest(a, state, bits1, 3);
        else
            acb_one(a);

        acb_poly_zeta_series(A, S, a, deflate, n1, bits2);
        acb_poly_zeta_series(B, S, a, deflate, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, FLINT_MIN(n1, n2));
        acb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(B, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* check zeta(s,a) = zeta(s,a+1) + a^(-s) */
        acb_poly_set_acb(D, a);
        acb_poly_log_series(D, D, n1, bits2);
        acb_poly_mullow(D, D, S, n1, bits2);
        acb_poly_neg(D, D);
        acb_poly_exp_series(D, D, n1, bits2);

        acb_add_ui(a1, a, 1, bits2);
        acb_poly_zeta_series(E, S, a1, deflate, n1, bits2);
        acb_poly_add(E, E, D, bits2);

        if (!acb_poly_overlaps(A, E))
        {
            flint_printf("FAIL (functional equation)\n\n");
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("E = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_zeta_series(S, S, a, deflate, n1, bits2);
        if (!acb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
        acb_poly_clear(D);
        acb_poly_clear(E);
        acb_poly_clear(F);
        acb_clear(a);
        acb_clear(a1);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

