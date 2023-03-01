/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("hardy_z_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        acb_poly_t S, A, B, C;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 8);
        n1 = 1 + n_randint(state, 8);
        n2 = 1 + n_randint(state, 8);

        do {
            q = 1 + n_randint(state, 15);
        } while (q % 4 == 2);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        do {
            dirichlet_char_index(chi, G, n_randint(state, G->phi_q));
        } while (!dirichlet_char_is_primitive(G, chi));

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);

        acb_poly_randtest(S, state, m, bits1, 3);

        acb_dirichlet_hardy_z_series(A, S, G, chi, n1, bits2);
        acb_poly_set(B, S);  /* aliasing */
        if (q == 1 && n_randint(state, 2))
            acb_poly_neg(B, B);
        acb_dirichlet_hardy_z_series(B, B, G, chi, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, FLINT_MIN(n1, n2));
        acb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(B, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

