/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "dirichlet.h"

/* index of the first difference, or -1 */
static slong
vec_diff(ulong * v, ulong * ref, ulong nv)
{
    ulong k;
    for (k = 0; k < nv; k++)
        if (ref[k] != v[k])
            return k;
    return -1;
}

static int
group_equal(const dirichlet_group_t G, const dirichlet_group_t H)
{
    slong k;

    if (G->q != H->q || G->q_even != H->q_even || G->rad_q != H->rad_q ||
        G->phi_q != H->phi_q || G->neven != H->neven || G->num != H->num ||
        G->expo != H->expo)
        return 0;

    for (k = 0; k < G->num; k++)
        if (G->P[k].p != H->P[k].p || G->P[k].e != H->P[k].e ||
            G->P[k].pe.n != H->P[k].pe.n || G->P[k].phi.n != H->P[k].phi.n ||
            G->generators[k] != H->generators[k] || G->PHI[k] != H->PHI[k])
            return 0;

    return 1;
}

TEST_FUNCTION_START(dirichlet_vec, state)
{
    ulong q;

    for (q = 2; q < 600; q ++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong * v1, * v2, nv;
        slong k;

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

            if ((k = vec_diff(v1, v2, nv)) >= 0 || v1[0] != DIRICHLET_CHI_NULL)
            {
                k = FLINT_MAX(k, 0);
                flint_printf("FAIL: chi_%wu(%wu,%wd) [mod %wu]\n", q, chi->n, k, G->expo);
                flint_printf("vec_loop      -> %wu\n", v1[k]);
                flint_printf("vec_primeloop -> %wu\n", v2[k]);
                flint_printf("pairing        = %wu\n", dirichlet_pairing(G, chi->n, k));
                flint_abort();
            }

            order = dirichlet_order_char(G, chi);
            dirichlet_chi_vec_loop_order(v1, G, chi, order, nv);
            dirichlet_chi_vec_primeloop_order(v2, G, chi, order, nv);

            if ((k = vec_diff(v1, v2, nv)) >= 0 || v1[0] != DIRICHLET_CHI_NULL)
            {
                k = FLINT_MAX(k, 0);
                flint_printf("FAIL: chi_%wu(%wu,%wd) [mod %wu]\n", q, chi->n, k, order);
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

    /* q = 1: the trivial character takes the value 1 (log 0) everywhere;
       nv = 1: only v[0] is written */
    for (q = 1; q < 8; q++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong v[4], w[4];
        slong i;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_one(chi, G);

        for (i = 0; i < 4; i++)
            v[i] = w[i] = 12345;

        dirichlet_chi_vec_loop(v, G, chi, 1);
        dirichlet_chi_vec_primeloop(w, G, chi, 1);

        if (v[0] != (q == 1 ? 0 : DIRICHLET_CHI_NULL) || v[0] != w[0] ||
            v[1] != 12345 || w[1] != 12345)
        {
            flint_printf("FAIL: nv = 1, q = %wu\n", q);
            flint_abort();
        }

        if (q == 1)
        {
            dirichlet_chi_vec(v, G, chi, 3);
            dirichlet_chi_vec_loop(w, G, chi, 3);

            if (v[0] != 0 || v[1] != 0 || v[2] != 0 ||
                w[0] != 0 || w[1] != 0 || w[2] != 0 || v[3] != 12345)
            {
                flint_printf("FAIL: q = 1\n");
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
    }

    /* subgroups agree with the groups of the same modulus */
    for (q = 1; q < 400; q++)
    {
        dirichlet_group_t G, H, K;
        ulong h;

        dirichlet_group_init(G, q);

        for (h = 1; h <= q; h++)
        {
            if (q % h)
                continue;

            /* garbage in all fields that dirichlet_subgroup_init sets */
            memset(H, 0xff, sizeof(dirichlet_group_t));
            dirichlet_subgroup_init(H, G, h);
            dirichlet_group_init(K, h);

            if (!group_equal(H, K))
            {
                flint_printf("FAIL: subgroup %wu of %wu\n", h, q);
                flint_abort();
            }

            dirichlet_group_clear(H);
            dirichlet_group_clear(K);
        }

        dirichlet_group_clear(G);
    }

#if FLINT_BITS == 64
    /* 16 components: q = 8 * 3 * 5 * ... * 47 */
    {
        dirichlet_group_t G, H, K;
        dirichlet_char_t chi;
        ulong v[60], w[60];
        slong i;

        q = 8;
        for (i = 3; i <= 47; i += 2)
            if (n_is_prime(i))
                q *= i;

        dirichlet_group_init(G, q);
        dirichlet_subgroup_init(H, G, q);
        dirichlet_group_init(K, q);

        if (G->num != 16 || !group_equal(H, K))
        {
            flint_printf("FAIL: 16 components\n");
            flint_abort();
        }

        /* the values on 1, ..., 59 are 0 or NULL for the trivial character */
        dirichlet_char_init(chi, G);
        dirichlet_char_one(chi, G);
        dirichlet_chi_vec(v, G, chi, 60);
        dirichlet_chi_vec_order(w, G, chi, 2, 60);

        for (i = 0; i < 60; i++)
        {
            ulong e = (n_gcd(i, q) == 1) ? 0 : DIRICHLET_CHI_NULL;

            if (v[i] != e || w[i] != e)
            {
                flint_printf("FAIL: 16 components, value at %wd\n", i);
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        dirichlet_group_clear(H);
        dirichlet_group_clear(K);
    }
#endif

    TEST_FUNCTION_END(state);
}
