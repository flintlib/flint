/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_gamma_lower, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a0, a1, b, z, w0, w1, t, u, enz;
        slong prec0, prec1;
        int regularized;

        acb_init(a0);
        acb_init(a1);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(t);
        acb_init(u);
        acb_init(enz);

        regularized = n_randint(state, 3);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_add_ui(a1, a0, 1, prec0);

        acb_hypgeom_gamma_lower(w0, a0, z, regularized, prec0);
        acb_hypgeom_gamma_lower(w1, a1, z, regularized, prec1);

        acb_neg(enz, z);
        acb_exp(enz, enz, prec0);

        /* recurrence relations */
        if (regularized == 2)
        {
            /* gamma^{*}(a,z) - exp(-z)/Gamma(a+1) - z gamma^{*}(a+1,z) = 0 */
            /* http://dlmf.nist.gov/8.8.E4 */
            acb_set(t, w0);
            acb_rgamma(u, a1, prec0);
            acb_submul(t, enz, u, prec0);
            acb_submul(t, z, w1, prec0);
        }
        else if (regularized == 1)
        {
            /* P(a,z) - exp(-z) z^a / Gamma(a+1) - P(a+1,z) = 0 */
            /* http://dlmf.nist.gov/8.8.E5 */
            acb_pow(u, z, a0, prec0);
            acb_rgamma(b, a1, prec0);
            acb_mul(u, u, b, prec0);
            acb_sub(t, w0, w1, prec0);
            acb_submul(t, enz, u, prec0);
        }
        else
        {
            /* a gamma(a,z) - exp(-z) z^a - gamma(a+1,z) = 0 */
            /* http://dlmf.nist.gov/8.8.E1 */
            acb_pow(u, z, a0, prec0);
            acb_mul(t, a0, w0, prec0);
            acb_submul(t, enz, u, prec0);
            acb_sub(t, t, w1, prec0);
        }

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: recurrence relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* identities relating lower and upper incomplete gamma functions */
        if (regularized == 0 || regularized == 1)
        {
            acb_t u0;
            acb_init(u0);
            acb_hypgeom_gamma_upper(u0, a0, z, regularized, prec0);

            acb_zero(t);

            if (regularized == 1)
            {
                /* P(s,z) + Q(s,z) - 1 = 0 */
                /* http://dlmf.nist.gov/8.2.E5 */
                acb_add(t, w0, u0, prec0);
                acb_sub_ui(t, t, 1, prec0);
            }
            else
            {
                /* gamma(s,z) + Gamma(s,z) - Gamma(s) = 0 */
                /* excludes non-positive integer values of s */
                /* http://dlmf.nist.gov/8.2.E3 */
                if (!acb_is_int(a0) || arb_is_positive(acb_realref(a0)))
                {
                    acb_gamma(b, a0, prec0);
                    acb_add(t, w0, u0, prec0);
                    acb_sub(t, t, b, prec0);
                }
            }

            if (!acb_contains_zero(t))
            {
                flint_printf("FAIL: lower plus upper\n\n");
                flint_printf("regularized = %d\n\n", regularized);
                flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
                flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
                flint_abort();
            }

            acb_clear(u0);
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(t);
        acb_clear(u);
        acb_clear(enz);
    }

    TEST_FUNCTION_END(state);
}
