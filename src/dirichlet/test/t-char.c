/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "dirichlet.h"

TEST_FUNCTION_START(dirichlet_char, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        dirichlet_group_t G;
        dirichlet_char_t x, y;
        ulong q, n, k, sum;
        slong ref;

        q = 1 + n_randint(state, 1000 * (1 + iter / 100));

        dirichlet_group_init(G, q);

        dirichlet_char_init(x, G);
        dirichlet_char_init(y, G);

        /* check group size and elements */
        dirichlet_char_one(x, G);
        sum = 1;

        for (n = 1; dirichlet_char_next(x, G) >= 0; n++)
            sum += x->n * x->n;

        if (FLINT_BITS == 64 || q < 1024)
        {
            /* use https://oeis.org/A053818 to check all elements
             * are gone through */
            ref = (q % 4 == 2) ? -2 : 1;
            for (k = (G->neven == 2); k < G->num; k++)
                ref = - ref * G->P[k].p;
            ref = ( G->phi_q * (2 * q * q + ref) ) / 6;

            if (n != G->phi_q)
            {
                flint_printf("FAIL: group size\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("phi(q) = %wu\n\n", G->phi_q);
                flint_printf("loop index = %wu\n\n", n);
                flint_abort();
            }
            if (sum != ref && q > 1)
            {
                flint_printf("FAIL: sum test\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("sum k^2 = %wu\n\n", ref);
                flint_printf("sum obtained = %wu\n\n", sum);
                flint_abort();
            }
        }

        if (q % 4 != 2)
        {
            dirichlet_char_first_primitive(x, G);
            for (n = 1; dirichlet_char_next_primitive(x, G) >= 0; n++);

            ref = dirichlet_group_num_primitive(G);
            if (n != ref)
            {
                flint_printf("FAIL: number of primitive elements\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("# primitive = %wu\n\n", ref);
                flint_printf("loop index = %wu\n\n", n);
                flint_abort();
            }

            /* some random elements, check log and exp */
            for (n = 0; n < 30; n++)
            {
                slong k;
                ulong m;

                for (m = 1; n_gcd(m, q) > 1; m = n_randint(state, q));
                dirichlet_char_log(x, G, m);

                if (m != _dirichlet_char_exp(x, G))
                {
                    flint_printf("FAIL: char log and exp\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("char = ");
                    dirichlet_char_print(G, x);
                    flint_printf("\n\nnumber = %wu\n\n", x->n);
                    flint_abort();
                }

                for (k = 0; k < G->num; k++)
                    x->log[k] = n_randint(state, G->P[k].phi.n);

                m = _dirichlet_char_exp(x, G);
                dirichlet_char_log(y, G, m);

                if (!dirichlet_char_eq_deep(G, x, y))
                {
                    flint_printf("FAIL: char exp and log\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("char = ");
                    dirichlet_char_print(G, x);
                    flint_printf("\n\nm = %wu\n\n", m);
                    flint_printf("log = ");
                    dirichlet_char_print(G, y);
                    flint_printf("\n\nnumber = %wu\n\n", y->n);
                    flint_abort();
                }

                dirichlet_char_next_primitive(x, G);
                m = x->n;

                if (m != _dirichlet_char_exp(x, G))
                {
                    flint_printf("FAIL: char number next primitive\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("char = ");
                    dirichlet_char_print(G, y);
                    flint_printf(", m = %wu\n\n", y->n);
                    flint_printf("next primitive = ");
                    dirichlet_char_print(G, x);
                    flint_printf(", m = %wu\n\n", m);
                    flint_printf("exp = %wu\n\n", x->n);
                    flint_abort();
                }
            }
        }

        dirichlet_char_clear(x);
        dirichlet_char_clear(y);
        dirichlet_group_clear(G);
    }

    TEST_FUNCTION_END(state);
}
