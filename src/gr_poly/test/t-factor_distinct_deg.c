/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_factor_distinct_deg, state)
{
    slong iter;

    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 3));
        gr_ctx_t ctx;
        gr_poly_t A, P, Q, G;
        gr_poly_vec_t fac, ed, irr;
        fmpz_vec_t degs;
        slong i, j, k, nfac, d, deg;
        int status = GR_SUCCESS;

        gr_ctx_init_random_finite_field(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(G, ctx);
        gr_poly_vec_init(fac, 0, ctx);
        gr_poly_vec_init(ed, 0, ctx);
        gr_poly_vec_init(irr, 0, ctx);
        fmpz_vec_init(degs, 0);

        /* build a product of distinct monic irreducible polynomials */
        nfac = 1 + n_randint(state, 5);
        status |= gr_poly_one(P, ctx);

        for (i = 0; i < nfac; i++)
        {
            deg = 1 + n_randint(state, 6);

            for (k = 0; k < 100; k++)
            {
                status |= gr_poly_randtest(A, state, deg, ctx);
                status |= gr_poly_set_coeff_ui(A, deg, 1, ctx);

                if (status == GR_SUCCESS && gr_poly_is_irreducible(A, ctx) == T_TRUE)
                    break;
            }

            if (k == 100)
                continue;

            /* skip repeated factors */
            for (j = 0; j < irr->length; j++)
                if (gr_poly_equal(A, irr->entries + j, ctx) != T_FALSE)
                    break;
            if (j < irr->length)
                continue;

            status |= gr_poly_vec_append(irr, A, ctx);
            status |= gr_poly_mul(P, P, A, ctx);
        }

        if (irr->length == 0)
            goto cleanup;

        /* random scalar multiple */
        {
            gr_ptr c;
            GR_TMP_INIT(c, ctx);
            status |= gr_randtest_not_zero(c, state, ctx);
            status |= gr_poly_mul_scalar(P, P, c, ctx);
            GR_TMP_CLEAR(c, ctx);
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (setup)\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        status = gr_poly_factor_distinct_deg(fac, degs, P, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (distinct_deg status)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_abort();
        }

        /* the product should be the monic version of P */
        status |= gr_poly_make_monic(Q, P, ctx);
        status |= gr_poly_one(G, ctx);
        for (i = 0; i < fac->length; i++)
            status |= gr_poly_mul(G, G, fac->entries + i, ctx);

        if (status != GR_SUCCESS || gr_poly_equal(G, Q, ctx) != T_TRUE)
        {
            flint_printf("FAIL (product)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("fac = %{gr_poly*}\n", fac->entries, fac->length, ctx);
            flint_printf("degs = %{fmpz*}\n", degs->entries, degs->length);
            flint_abort();
        }

        /* each factor should split into irreducibles of the given degree */
        for (i = 0; i < fac->length; i++)
        {
            d = fmpz_get_si(degs->entries + i);
            deg = fac->entries[i].length - 1;

            if (d <= 0 || deg % d != 0)
            {
                flint_printf("FAIL (degree)\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("fac = %{gr_poly*}\n", fac->entries, fac->length, ctx);
                flint_printf("degs = %{fmpz*}\n", degs->entries, degs->length);
                flint_abort();
            }

            /* distinct degrees */
            for (j = 0; j < i; j++)
            {
                if (fmpz_equal(degs->entries + i, degs->entries + j))
                {
                    flint_printf("FAIL (repeated degree)\n");
                    gr_ctx_println(ctx);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("fac = %{gr_poly*}\n", fac->entries, fac->length, ctx);
                    flint_printf("degs = %{fmpz*}\n", degs->entries, degs->length);
                    flint_abort();
                }
            }

            status = gr_poly_factor_equal_deg(ed, fac->entries + i, d, ctx);

            status |= gr_poly_one(G, ctx);
            for (j = 0; j < ed->length; j++)
            {
                status |= gr_poly_mul(G, G, ed->entries + j, ctx);

                if (ed->entries[j].length - 1 != d)
                {
                    flint_printf("FAIL (equal_deg degree)\n");
                    gr_ctx_println(ctx);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("f = "); gr_poly_print(fac->entries + i, ctx); flint_printf("\n");
                    flint_printf("ed = %{gr_poly*}\n", ed->entries, ed->length, ctx);
                    flint_abort();
                }

                for (k = 0; k < irr->length; k++)
                    if (gr_poly_equal(ed->entries + j, irr->entries + k, ctx) == T_TRUE)
                        break;

                if (k == irr->length)
                {
                    flint_printf("FAIL (equal_deg factor not in list)\n");
                    gr_ctx_println(ctx);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("f = "); gr_poly_print(fac->entries + i, ctx); flint_printf("\n");
                    flint_printf("ed = %{gr_poly*}\n", ed->entries, ed->length, ctx);
                    flint_printf("irr = %{gr_poly*}\n", irr->entries, irr->length, ctx);
                    flint_abort();
                }
            }

            if (status != GR_SUCCESS || ed->length != deg / d ||
                gr_poly_equal(G, fac->entries + i, ctx) != T_TRUE)
            {
                flint_printf("FAIL (equal_deg product)\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("f = "); gr_poly_print(fac->entries + i, ctx); flint_printf("\n");
                flint_printf("ed = %{gr_poly*}\n", ed->entries, ed->length, ctx);
                flint_abort();
            }

            /* test the probabilistic splitting */
            if (deg > d)
            {
                status = gr_poly_factor_equal_deg_prob(G, state, fac->entries + i, d, ctx);

                if (status != GR_SUCCESS ||
                    !(G->length == 1 || (G->length > 1 && G->length < fac->entries[i].length &&
                        (G->length - 1) % d == 0)))
                {
                    flint_printf("FAIL (equal_deg_prob)\n");
                    gr_ctx_println(ctx);
                    flint_printf("f = "); gr_poly_print(fac->entries + i, ctx); flint_printf("\n");
                    flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                    flint_abort();
                }

                if (G->length > 1)
                {
                    status = gr_poly_divrem(A, Q, fac->entries + i, G, ctx);
                    if (status != GR_SUCCESS || Q->length != 0)
                    {
                        flint_printf("FAIL (equal_deg_prob divisibility)\n");
                        gr_ctx_println(ctx);
                        flint_printf("f = "); gr_poly_print(fac->entries + i, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }

cleanup:
        gr_poly_clear(A, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_vec_clear(fac, ctx);
        gr_poly_vec_clear(ed, ctx);
        gr_poly_vec_clear(irr, ctx);
        fmpz_vec_clear(degs);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
