/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* adapted from code generated using Claude Fable 5 */

#include "test_helpers.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

TEST_GR_FUNCTION_START(gr_ore_poly_indicial_polynomial, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t Pol;
        gr_ore_poly_ctx_t Ore;
        gr_ore_poly_t op;
        gr_poly_struct * ind, * pol;
        int status = GR_SUCCESS;
        int differential = 0;

        if (!n_randint(state, 8))
            gr_ore_poly_ctx_init_randtest2(Pol, Ore, state);
        else
        {
            gr_ctx_init_random_poly(Pol, state);
            ore_algebra_t alg = n_randint(state, 2)
                ? ORE_ALGEBRA_DERIVATIVE
                : ORE_ALGEBRA_EULER_DERIVATIVE;
            gr_ore_poly_ctx_init(Ore, Pol, 0, alg);
            differential = 1;
        }

        gr_ore_poly_init(op, Ore);
        ind = (gr_poly_struct *) gr_heap_init(Pol);
        pol = (gr_poly_struct *) gr_heap_init(Pol);

        status |= gr_ore_poly_randtest(op, state, 1 + n_randint(state, 10), Ore);

        status |= gr_ore_poly_indicial_polynomial(ind, op, Ore);

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        if (status == GR_SUCCESS && differential)
        {
            gr_ctx_struct * Cst;
            status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, Ore);
            GR_MUST_SUCCEED(status);

            gr_ptr c = gr_heap_init(Cst);
            gr_ptr nn = gr_heap_init(Cst);

            for (slong k = 0; k < 10; k++)
            {
                /* check that op(x^n) = ind(n)*x^{n+v} + ... */
                slong n = n_randint(state, 1000);
                status |= gr_gen(pol, Pol);
                status |= gr_pow_ui(pol, pol, n, Pol);
                status |= gr_ore_poly_apply(pol, op, pol, Ore);

                status |= gr_set_ui(nn, n, Cst);
                status |= gr_poly_evaluate(c, ind, nn, Cst);

                for (slong i = 0; i < pol->length; i++)
                {
                    gr_ptr c1 = gr_poly_coeff_ptr(pol, i, Cst);
                    if (gr_equal(c1, c, Cst) != T_FALSE)
                        break;
                    else if (gr_is_zero(c1, Cst) != T_FALSE)
                        continue;

                    flint_printf("FAIL: valuation coefficient\n\n");
                    flint_printf("Ore = %{gr_ctx}\n", Ore);
                    flint_printf("op = %{gr}\n", op, Ore);
                    flint_printf("n = %wd\n", n);
                    flint_abort();
                }
            }

            gr_heap_clear(c, Cst);
            gr_heap_clear(nn, Cst);
            gr_ctx_clear(Cst);
        }

        gr_clear(ind, Pol);
        gr_ore_poly_clear(op, Ore);
        gr_ore_poly_ctx_clear(Ore);
        gr_ctx_clear(Pol);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
