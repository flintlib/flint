/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_roots_finite_field, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, P, Q, R;
        gr_vec_t roots, roots2;
        fmpz_vec_t mult, mult2;
        gr_ptr c, t;
        slong i, j, nroots, deg;
        ulong e;
        int status = GR_SUCCESS;

        gr_ctx_init_random_finite_field(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(R, ctx);
        gr_vec_init(roots, 0, ctx);
        gr_vec_init(roots2, 0, ctx);
        fmpz_vec_init(mult, 0);
        fmpz_vec_init(mult2, 0);
        c = gr_heap_init(ctx);
        t = gr_heap_init(ctx);

        /* P = random poly times product of (x - r)^e */
        status |= gr_poly_randtest(P, state, 1 + n_randint(state, 8), ctx);
        if (P->length == 0)
            status |= gr_poly_one(P, ctx);

        nroots = n_randint(state, 5);
        for (i = 0; i < nroots; i++)
        {
            status |= gr_randtest(c, state, ctx);
            status |= gr_neg(c, c, ctx);
            status |= gr_poly_gen(A, ctx);
            status |= gr_poly_set_coeff_scalar(A, 0, c, ctx);
            status |= gr_poly_pow_ui(A, A, 1 + n_randint(state, 3), ctx);
            status |= gr_poly_mul(P, P, A, ctx);
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (setup)\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        status = gr_poly_roots_finite_field(roots, mult, P, 0, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (status)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_abort();
        }

        /* Q = P / prod (x - r_i)^m_i should have no roots among the r_i
           and the r_i must be distinct */
        status |= gr_poly_set(Q, P, ctx);
        for (i = 0; i < roots->length; i++)
        {
            e = fmpz_get_ui(mult->entries + i);
            status |= gr_poly_gen(A, ctx);
            status |= gr_neg(c, gr_vec_entry_srcptr(roots, i, ctx), ctx);
            status |= gr_poly_set_coeff_scalar(A, 0, c, ctx);
            status |= gr_poly_pow_ui(A, A, e, ctx);
            status |= gr_poly_divrem(Q, R, Q, A, ctx);

            if (status != GR_SUCCESS || R->length != 0)
            {
                flint_printf("FAIL (divisibility)\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("roots = %{gr*}\n", roots->entries, roots->length, ctx);
                flint_printf("mult = %{fmpz*}\n", mult->entries, mult->length);
                flint_abort();
            }

            for (j = 0; j < i; j++)
            {
                if (gr_equal(gr_vec_entry_srcptr(roots, i, ctx), gr_vec_entry_srcptr(roots, j, ctx), ctx) != T_FALSE)
                {
                    flint_printf("FAIL (repeated root)\n");
                    gr_ctx_println(ctx);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("roots = %{gr*}\n", roots->entries, roots->length, ctx);
                    flint_abort();
                }
            }
        }

        for (i = 0; i < roots->length; i++)
        {
            status |= gr_poly_evaluate(t, Q, gr_vec_entry_srcptr(roots, i, ctx), ctx);

            if (status != GR_SUCCESS || gr_is_zero(t, ctx) != T_FALSE)
            {
                flint_printf("FAIL (multiplicity)\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("roots = %{gr*}\n", roots->entries, roots->length, ctx);
                flint_printf("mult = %{fmpz*}\n", mult->entries, mult->length);
                flint_abort();
            }
        }

        /* compare with the ring's own root finder */
        status = gr_poly_roots(roots2, mult2, P, 0, ctx);

        if (status == GR_SUCCESS && roots2->length != roots->length)
        {
            flint_printf("FAIL (number of roots)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("roots = %{gr*}\n", roots->entries, roots->length, ctx);
            flint_printf("roots2 = %{gr*}\n", roots2->entries, roots2->length, ctx);
            flint_abort();
        }

        /* check the deflation routines */
        deg = gr_poly_deflation(P, ctx);
        status = gr_poly_deflate(A, P, deg, ctx);
        status |= gr_poly_inflate(Q, A, deg, ctx);

        if (status != GR_SUCCESS || deg == 0 || gr_poly_equal(P, Q, ctx) != T_TRUE
                || (P->length > 1 && ((P->length - 1) % deg != 0 || A->length - 1 != (P->length - 1) / deg)))
        {
            flint_printf("FAIL (deflation)\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
            flint_printf("deflation = %wd\n", deg);
            flint_abort();
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(R, ctx);
        gr_vec_clear(roots, ctx);
        gr_vec_clear(roots2, ctx);
        fmpz_vec_clear(mult);
        fmpz_vec_clear(mult2);
        gr_heap_clear(c, ctx);
        gr_heap_clear(t, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
