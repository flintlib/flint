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

static void
check_eq(ulong p1, ulong p2, ulong q, ulong m, char * fct1, char * fct2)
{
    if (p1 != p2)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("chi_%wu(%wu,.)\n\n", q, m);
        flint_printf("%s = %wu\n\n", fct1, p1);
        flint_printf("%s = %wu\n\n", fct2, p2);
        flint_abort();
    }
}

static ulong
random_divisor(flint_rand_t state, const dirichlet_group_t G)
{
    ulong d;
    slong k;
    d = 1;
    for (k = (G->neven == 2); k < G->num; k++)
        d *= n_pow(G->P[k].p, n_randint(state, G->P[k].e));
    return d;
}

TEST_FUNCTION_START(dirichlet_properties, state)
{
    slong iter, bits;

    for (bits = 5; bits <= 30; bits += 5)
    {

        for (iter = 0; iter < 50; iter++)
        {
            dirichlet_group_t G;
            dirichlet_char_t chi, psi;
            ulong q, iter2;

            q = 2 + n_randint(state, 1 << bits);

            dirichlet_group_init(G, q);
            dirichlet_char_init(chi, G);
            dirichlet_char_init(psi, G);

            /* check number char properties */
            for (iter2 = 0; iter2 < 100; iter2++)
            {
                ulong m, n;
                ulong p1, p2, pairing, cm, cn, q2, q3;
                dirichlet_group_t G2, G3;
                dirichlet_char_t chi2, chi3;

                if (iter2 == 50)
                    dirichlet_group_dlog_precompute(G, 5);

                /* one random character */
                do
                    m = n_randint(state, q);
                while (n_gcd(q, m) > 1);

                dirichlet_char_log(chi, G, m);

                p1 = dirichlet_order_ui(G, m);
                p2 = dirichlet_order_char(G, chi);
                check_eq(p1, p2, q, m, "order m", "order chi");

                p1 = dirichlet_conductor_ui(G, m);
                p2 = dirichlet_conductor_char(G, chi);
                check_eq(p1, p2, q, m, "conductor m", "conductor chi");

                p1 = dirichlet_parity_ui(G, m);
                p2 = dirichlet_parity_char(G, chi);
                check_eq(p1, p2, q, m, "parity m", "parity chi");

                p1 = dirichlet_char_is_real(G, chi);
                p2 = (dirichlet_order_char(G, chi) <= 2);
                check_eq(p1, p2, q, m, "is_real", "(order <= 2)");

                /* check index */
                p1 = dirichlet_index_char(G, chi);
                dirichlet_char_index(psi, G, p1);

                if (!dirichlet_char_eq_deep(G, chi, psi))
                {
                    flint_printf("FAIL: index\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("chi = "); dirichlet_char_print(G, chi);
                    flint_printf("\n\nindex(chi) = %wu\n\n", p1);
                    flint_printf("psi(index) = %wu\n\n", psi->n);
                    flint_printf("psi = "); dirichlet_char_print(G, psi);
                    flint_printf("\n\n");
                    flint_abort();
                }

                /* lift to higher modulus */
                do {
                    q2 = q * (1 + n_randint(state, 100));
                } while (q2 % q); /* in case of overflow */

                dirichlet_group_init(G2, q2);
                dirichlet_char_init(chi2, G2);
                dirichlet_char_lift(chi2, G2, chi, G);

                p1 = dirichlet_conductor_char(G, chi);
                p2 = dirichlet_conductor_char(G2, chi2);
                check_eq(p1, p2, q, m, "conductor chi", "conductor lift");

                p1 = dirichlet_order_char(G, chi);
                p2 = dirichlet_order_char(G2, chi2);
                check_eq(p1, p2, q, m, "order chi", "order lift");

                /* and lower */

                dirichlet_char_lower(psi, G, chi2, G2);
                if (!dirichlet_char_eq_deep(G, chi, psi))
                {
                    flint_printf("FAIL: lift and lower back\n\n");
                    flint_printf("q = %wu\n\nchi = ", q);
                    dirichlet_char_print(G, chi);
                    flint_printf("\n\nq2 = %wu\n\nchi2 = ", q2);
                    dirichlet_char_print(G2, chi2);
                    flint_printf("\n\nq = %wu\n\npsi = ", q);
                    dirichlet_char_print(G, psi);
                    flint_printf("\n\n");
                    flint_abort();
                }

                /* choose q3 s.t. conductor | q3 | q */
                q3 = dirichlet_conductor_char(G, chi) * random_divisor(state, G);
                q3 = n_gcd(q, q3);

                /* discard if overflow */
                if (q3 % p1 == 0)
                {
                    dirichlet_group_init(G3, q3);
                    dirichlet_char_init(chi3, G3);
                    dirichlet_char_lower(chi3, G3, chi2, G2);

                    p1 = dirichlet_conductor_char(G, chi);
                    p2 = dirichlet_conductor_char(G3, chi3);
                    check_eq(p1, p2, q, m, "conductor chi", "conductor lower");

                    p1 = dirichlet_order_char(G, chi);
                    p2 = dirichlet_order_char(G3, chi3);
                    check_eq(p1, p2, q, m, "order chi", "order lower");

                    dirichlet_char_clear(chi3);
                    dirichlet_group_clear(G3);
                }

                dirichlet_char_clear(chi2);
                dirichlet_group_clear(G2);

                /* another random character */
                do
                    n = n_randint(state, q);
                while (n_gcd(q, n) > 1);

                dirichlet_char_log(psi, G, n);
                pairing = dirichlet_pairing(G, m, n);
                cn = dirichlet_chi(G, chi, n);
                cm = dirichlet_chi(G, psi, m);

                if (pairing != cn || pairing != cm)
                {
                    flint_printf("FAIL: pairing\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("n = %wu\n\n", n);
                    flint_printf("chi(m,n) = %wu\n\n", pairing);
                    flint_printf("chi(m)(n) = %wu\n\n", cn);
                    flint_printf("chi(n)(m) = %wu\n\n", cm);
                    flint_abort();
                }

            }

            dirichlet_group_dlog_clear(G);

            dirichlet_char_clear(chi);
            dirichlet_char_clear(psi);
            dirichlet_group_clear(G);
        }
    }

    TEST_FUNCTION_END(state);
}
