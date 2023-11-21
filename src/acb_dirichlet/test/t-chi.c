/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_chi, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t zn1, zn2, zn1n2, zn1zn2;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, m, n1, n2, iter2;
        int res;

        q = 1 + n_randint(state, 1000);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        acb_init(zn1);
        acb_init(zn2);
        acb_init(zn1n2);
        acb_init(zn1zn2);

        /* check chi(n1) chi(n2) = chi(n1 n2) */
        for (iter2 = 0; iter2 < 10; iter2++)
        {
            do {
                m = 1 + n_randint(state, q);
            } while (n_gcd(q, m) != 1);

            dirichlet_char_log(chi, G, m);

            n1 = n_randint(state, 1000);
            n2 = n_randint(state, 1000);

            acb_dirichlet_chi(zn1, G, chi, n1, 53);
            acb_dirichlet_pairing(zn2, G, m, n1, 53);
            if (!acb_overlaps(zn1, zn2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("n = %wu\n\n", n1);
                flint_printf("char = "); acb_printd(zn1, 15); flint_printf("\n\n");
                flint_printf("pairing = "); acb_printd(zn2, 15); flint_printf("\n\n");
                dirichlet_char_print(G, chi);
                dirichlet_char_log(chi, G, m);
                flint_printf("log(m) = "); dirichlet_char_print(G, chi);
                dirichlet_char_log(chi, G, n1);
                flint_printf("log(n1) = "); dirichlet_char_print(G, chi);
                flint_abort();
            }

            acb_dirichlet_pairing(zn2, G, m, n2, 53);
            acb_dirichlet_pairing(zn1n2, G, m, n1 * n2, 53);
            acb_mul(zn1zn2, zn1, zn2, 53);

            if (!acb_overlaps(zn1n2, zn1zn2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("n1 = %wu\n\n", n1);
                flint_printf("n2 = %wu\n\n", n2);
                flint_printf("zn1 = "); acb_printd(zn1, 15); flint_printf("\n\n");
                flint_printf("zn2 = "); acb_printd(zn2, 15); flint_printf("\n\n");
                flint_printf("zn1n2 = "); acb_printd(zn1n2, 15); flint_printf("\n\n");
                flint_printf("zn1zn2 = "); acb_printd(zn1zn2, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (iter % 10 == 0)
        {
            /* check orthogonality */
            acb_zero(zn1);
            n1 = n_randint(state, 1000);
            for (m = 1; m <= q; m++)
            {
                if (n_gcd(q, m) == 1)
                {
                    acb_dirichlet_pairing(zn2, G, m, n1, 53);
                    acb_add(zn1, zn1, zn2, 53);
                }
            }

            if (n1 % q == 1 % q)
                res = arb_contains_si(acb_realref(zn1), n_euler_phi(q)) &&
                    arb_contains_zero(acb_imagref(zn1));
            else
                res = acb_contains_zero(zn1);

            if (!res)
            {
                flint_printf("FAIL: orthogonality\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("phi = %wu\n\n", n_euler_phi(q));
                flint_printf("n1 = %wu\n\n", n1);
                flint_printf("zn1 = "); acb_printd(zn1, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
        acb_clear(zn1);
        acb_clear(zn2);
        acb_clear(zn1n2);
        acb_clear(zn1zn2);
    }

    TEST_FUNCTION_END(state);
}
