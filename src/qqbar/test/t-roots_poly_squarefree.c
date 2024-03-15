/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_poly.h"

TEST_FUNCTION_START(qqbar_roots_poly_squarefree, state)
{
    slong iter;
    slong count_success = 0;

    slong deg_limit = 200;
    slong bits_limit = 1000;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_poly_t poly, fac;
        gr_ctx_t ctx;
        qqbar_ptr roots, computed_roots;
        int success, status;
        qqbar_t c;
        slong deg, i, j;

        deg = n_randint(state, 4);

        roots = _qqbar_vec_init(deg);
        computed_roots = _qqbar_vec_init(deg);
        qqbar_init(c);

        gr_ctx_init_complex_qqbar(ctx);
        _gr_ctx_qqbar_set_limits(ctx, deg_limit, bits_limit);

        gr_poly_init(poly, ctx);
        gr_poly_init(fac, ctx);

        for (i = 0; i < deg; i++)
        {
            if (n_randint(state, 10) == 0)
                qqbar_randtest(roots + i, state, 4, 4);
            else
                qqbar_randtest(roots + i, state, 2, 4);

            for (j = 0; j < i; j++)
            {
                if (qqbar_equal(roots + i, roots + j))
                {
                    i--;
                    break;
                }
            }
        }

        if (n_randint(state, 5))
        {
            qqbar_one(c);
        }
        else
        {
            do {
                qqbar_randtest(c, state, 2, 4);
            } while (qqbar_is_zero(c));
        }

        status = GR_SUCCESS;
        status |= gr_poly_set_scalar(poly, c, ctx);

        for (i = 0; i < deg; i++)
        {
            status |= gr_poly_set_scalar(fac, roots + i, ctx);
            status |= gr_poly_neg(fac, fac, ctx);
            status |= gr_poly_set_coeff_si(fac, 1, 1, ctx);
            status |= gr_poly_mul(poly, poly, fac, ctx);
        }

        if (status == GR_SUCCESS)
        {
            success = _qqbar_roots_poly_squarefree(computed_roots, poly->coeffs, poly->length, deg_limit, bits_limit);
            count_success += success;

            if (success)
            {
                for (i = 0; i < deg; i++)
                {
                    int found = 0;

                    for (j = 0; j < deg; j++)
                    {
                        if (qqbar_equal(roots + i, computed_roots + j))
                        {
                            found = 1;
                            break;
                        }
                    }

                    if (!found)
                    {
                        flint_printf("FAIL\n\n");
                        flint_printf("deg = %wd\n", deg);
                        flint_printf("roots:\n");

                        for (j = 0; j < deg; j++)
                        {
                            qqbar_print(roots + i);
                            flint_printf("\n");
                        }

                        flint_printf("computed roots:\n");
                        for (j = 0; j < deg; j++)
                        {
                            qqbar_print(roots + i);
                            flint_printf("\n");
                        }

                        flint_abort();
                    }
                }
            }

/*            flint_printf("total: %wd / %wd\n", count_success, iter);  */
        }

        gr_poly_clear(poly, ctx);
        gr_poly_clear(fac, ctx);
        gr_ctx_clear(ctx);
        _qqbar_vec_clear(roots, deg);
        _qqbar_vec_clear(computed_roots, deg);
        qqbar_clear(c);
    }

    if (iter > 100 && count_success < 0.1 * iter)
    {
        flint_printf("FAIL: only %wd / %wd succeeded\n", count_success, iter);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
