/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_thetanull, state)
{
    slong prec = 64;
    ulong q;

    /* check the only theta functions
     *   theta(chi) = sum chi(k)* k^odd * exp(-Pi * k^2 / q)
     * vanishing at 1 correspond to two specific
     * characters of moduli 300 and 600 + conjugates
     */

    for (q = 3; q < 1000; q ++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong * v, nv, k;

        acb_t sum;
        acb_ptr z;

        arb_t eq;
        arb_ptr t, kt, tt;

        if (q % 4 == 2)
            /* no primitive character mod q */
            continue;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);

        z = _acb_vec_init(G->expo);
        _acb_vec_unit_roots(z, G->expo, G->expo, prec);

        nv = acb_dirichlet_theta_length_d(q, 1, prec);
        v = flint_malloc(nv * sizeof(ulong));

        arb_init(eq);
        arb_const_pi(eq, prec);
        arb_div_ui(eq, eq, q, prec);
        arb_neg(eq, eq);
        arb_exp(eq, eq, prec);

        t = _arb_vec_init(nv);
        acb_dirichlet_arb_quadratic_powers(t, nv, eq, prec);

        kt = _arb_vec_init(nv);
        for (k = 1; k < nv; k++)
            arb_mul_ui(kt + k, t + k, k, prec);

        /* theta function on primitive characters */
        acb_init(sum);
        dirichlet_char_first_primitive(chi, G);

        do {

            acb_zero(sum);
            dirichlet_chi_vec(v, G, chi, nv);

            tt = dirichlet_parity_char(G, chi) ? kt : t;

            for (k = 1; k < nv; k++)
                if (v[k] != DIRICHLET_CHI_NULL)
                    acb_addmul_arb(sum, z + v[k], tt + k, prec);

            if ((q == 300 && (chi->n == 71 || chi->n == 131))
                    || (q == 600 && (chi->n == 11 || chi->n == 491)))
            {
                if (!acb_contains_zero(sum))
                {
                flint_printf("FAIL: Theta(chi_%wu(%wu))=", q, chi->n);
                acb_printd(sum, 10);
                flint_printf("\n");
                dirichlet_char_print(G, chi);
                flint_printf("\n");
                flint_abort();
                }
            }
            else if (acb_contains_zero(sum))
            {
                flint_printf("FAIL: Theta(chi_%wu(%wu))=", q, chi->n);
                acb_printd(sum, 10);
                flint_printf("\n");
                dirichlet_char_print(G, chi);
                flint_printf("\n");
                flint_abort();
            }

        } while (dirichlet_char_next_primitive(chi, G) >= 0);

        _acb_vec_clear(z, G->expo);
        _arb_vec_clear(kt, nv);
        _arb_vec_clear(t, nv);
        acb_clear(sum);
        arb_clear(eq);
        flint_free(v);
        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
    }

    TEST_FUNCTION_END(state);
}
