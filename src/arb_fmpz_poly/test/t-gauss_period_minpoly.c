/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_fmpz_poly.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(arb_fmpz_poly_gauss_period_minpoly, state)
{
    slong iter;

    {
        slong prec;
        ulong n, q;
        fmpz_poly_t pol;
        fmpz_poly_init(pol);

        /* test q = 0 separately; see issue #194 */
        for (n = 0; n < 100; n++)
        {
            fmpz_poly_one(pol);
            arb_fmpz_poly_gauss_period_minpoly(pol, 0, n);
            if (!fmpz_poly_is_zero(pol))
            {
                flint_printf("FAIL (q = 0)\n");
                flint_abort();
            }
        }

        for (q = 1; q < 1000; q++)
        {
            acb_dirichlet_roots_t zeta;
            prec = 100 + n_randint(state, 500);

            if (n_is_prime(q))
                acb_dirichlet_roots_init(zeta, q, 30, prec);

            for (n = 0; n < 1000; n++)
            {
                arb_fmpz_poly_gauss_period_minpoly(pol, q, n);

                if (!fmpz_poly_is_zero(pol))
                {
                    ulong k, g, gk, e, d;
                    acb_t t, u;

                    acb_init(t);
                    acb_init(u);
                    d = (q - 1) / n;
                    g = n_primitive_root_prime(q);

                    for (iter = 0; iter < 3; iter++)
                    {
                        k = n_randint(state, n);
                        gk = n_powmod2(g, k, 2 * q);
                        acb_zero(u);

                        for (e = 0; e < d; e++)
                        {
                            acb_dirichlet_root(t, zeta, n_mulmod2_preinv(gk,
                                n_powmod2(g, n * e, 2 * q), 2 * q, n_preinvert_limb(2 * q)), prec);
                            acb_add(u, u, t, prec);
                        }

                        arb_fmpz_poly_evaluate_acb(t, pol, u, prec);

                        if (!acb_contains_zero(t) || fmpz_poly_degree(pol) != n)
                        {
                            flint_printf("FAIL\n");
                            flint_printf("q = %wu, n = %wu, k = %wu\n\n", q, n, k);
                            fmpz_poly_print(pol); flint_printf("\n\n");
                            acb_printn(u, 30, 0); flint_printf("\n\n");
                            acb_printn(t, 30, 0); flint_printf("\n\n");
                            flint_abort();
                        }
                    }

                    acb_clear(t);
                    acb_clear(u);
                }
            }

            if (n_is_prime(q))
                acb_dirichlet_roots_clear(zeta);
        }

        fmpz_poly_clear(pol);
    }

    TEST_FUNCTION_END(state);
}
