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

static ulong
vec_diff(ulong * v, ulong * ref, ulong nv)
{
    ulong k;
    for (k = 1; k < nv; k++)
        if (ref[k] != v[k])
            return k;
    return 0;
}

TEST_FUNCTION_START(dirichlet_vec, state)
{
    ulong q;

    for (q = 2; q < 600; q ++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong * v1, * v2, nv, k;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);

        nv = 100;
        v1 = flint_malloc(nv * sizeof(ulong));
        v2 = flint_malloc(nv * sizeof(ulong));

        dirichlet_char_one(chi, G);

        do {

            ulong order;

            dirichlet_chi_vec_loop(v1, G, chi, nv);
            dirichlet_chi_vec_primeloop(v2, G, chi, nv);

            if ((k = vec_diff(v1, v2, nv)))
            {
                flint_printf("FAIL: chi_%wu(%wu,%wu) [mod %wu]\n", q, chi->n, k, G->expo);
                flint_printf("vec_loop      -> %wu\n", v1[k]);
                flint_printf("vec_primeloop -> %wu\n", v2[k]);
                flint_printf("pairing        = %wu\n", dirichlet_pairing(G, chi->n, k));
                flint_abort();
            }

            order = dirichlet_order_char(G, chi);
            dirichlet_chi_vec_loop_order(v1, G, chi, order, nv);
            dirichlet_chi_vec_primeloop_order(v2, G, chi, order, nv);

            if ((k = vec_diff(v1, v2, nv)))
            {
                flint_printf("FAIL: chi_%wu(%wu,%wu) [mod %wu]\n", q, chi->n, k, order);
                flint_printf("vec_loop      -> %wu\n", v1[k]);
                flint_printf("vec_primeloop -> %wu\n", v2[k]);
                flint_printf("pairing        = %wu mod %wu\n", dirichlet_pairing(G, chi->n, k), G->expo);
                flint_abort();
            }

        } while (dirichlet_char_next(chi, G) >= 0);

        flint_free(v1);
        flint_free(v2);
        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
    }

    TEST_FUNCTION_END(state);
}
