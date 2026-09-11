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

/* Random finite field; with some probability, a small one which gives
   many repeated factors. */
static void
_init_ctx(gr_ctx_t ctx, flint_rand_t state)
{
    gr_ctx_init_random_finite_field(ctx, state);
}

TEST_FUNCTION_START(gr_poly_factor_finite_field, state)
{
    slong iter;

    for (iter = 0; iter < 150 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 3));
        gr_ctx_t ctx;
        gr_poly_t A, P, Q;
        gr_ptr c;
        gr_poly_vec_t fac;
        fmpz_vec_t exp;
        slong i, j, nfac, len_bound, exp_bound, deg, deg_sum;
        int algorithm;
        int status = GR_SUCCESS;

        _init_ctx(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_vec_init(fac, 0, ctx);
        fmpz_vec_init(exp, 0);
        c = gr_heap_init(ctx);

        nfac = n_randint(state, 5);
        len_bound = 1 + n_randint(state, 12);
        exp_bound = 1 + n_randint(state, 4);

        /* random product of random polynomials with random exponents */
        status |= gr_poly_one(P, ctx);
        for (i = 0; i < nfac; i++)
        {
            do
            {
                status |= gr_poly_randtest(A, state, len_bound, ctx);
            }
            while (A->length == 0);

            /* occasionally, build x^m - a to stress deflation and
               repeated factors */
            if (n_randint(state, 4) == 0)
            {
                status |= gr_poly_zero(A, ctx);
                status |= gr_poly_set_coeff_ui(A, 1 + n_randint(state, 12), 1, ctx);
                status |= gr_randtest(c, state, ctx);
                status |= gr_poly_set_coeff_scalar(A, 0, c, ctx);
            }

            status |= gr_poly_pow_ui(A, A, 1 + n_randint(state, exp_bound), ctx);
            status |= gr_poly_mul(P, P, A, ctx);
        }

        if (n_randint(state, 4) == 0)
            status |= gr_poly_randtest(P, state, 1 + n_randint(state, 40), ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (setup)\n");
            gr_ctx_println(ctx);
            flint_abort();
        }

        algorithm = n_randint(state, 4);

        status = _gr_poly_factor_finite_field(c, fac, exp, P, algorithm, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL (unexpected status %d)\n\n", status);
            gr_ctx_println(ctx);
            flint_printf("algorithm = %d\n", algorithm);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_abort();
        }

        /* check product */
        deg_sum = 0;
        status |= gr_poly_set_scalar(Q, c, ctx);
        for (i = 0; i < fac->length; i++)
        {
            deg = fac->entries[i].length - 1;

            if (!fmpz_fits_si(exp->entries + i) || fmpz_sgn(exp->entries + i) <= 0)
            {
                flint_printf("FAIL (bad exponent)\n\n");
                flint_abort();
            }

            deg_sum += deg * fmpz_get_si(exp->entries + i);

            status |= gr_poly_pow_ui(A, fac->entries + i, fmpz_get_ui(exp->entries + i), ctx);
            status |= gr_poly_mul(Q, Q, A, ctx);
        }

        if (status != GR_SUCCESS || gr_poly_equal(P, Q, ctx) != T_TRUE ||
            (P->length >= 1 && deg_sum != P->length - 1))
        {
            flint_printf("FAIL (product)\n\n");
            gr_ctx_println(ctx);
            flint_printf("algorithm = %d\n", algorithm);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
            flint_printf("c = %{gr}\n\n", c, ctx);
            flint_printf("fac = %{gr_poly*}\n\n", fac->entries, fac->length, ctx);
            flint_printf("exp = %{fmpz*}\n\n", exp->entries, exp->length);
            flint_abort();
        }

        /* check that the factors are monic, irreducible and distinct */
        for (i = 0; i < fac->length; i++)
        {
            gr_poly_struct * f = fac->entries + i;

            if (f->length < 2 || gr_is_one(gr_poly_coeff_srcptr(f, f->length - 1, ctx), ctx) != T_TRUE)
            {
                flint_printf("FAIL (not monic)\n\n");
                gr_ctx_println(ctx);
                flint_printf("algorithm = %d\n", algorithm);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("fac = %{gr_poly*}\n\n", fac->entries, fac->length, ctx);
                flint_abort();
            }

            if (gr_poly_is_irreducible(f, ctx) != T_TRUE)
            {
                flint_printf("FAIL (not irreducible)\n\n");
                gr_ctx_println(ctx);
                flint_printf("algorithm = %d\n", algorithm);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
                flint_printf("fac = %{gr_poly*}\n\n", fac->entries, fac->length, ctx);
                flint_abort();
            }

            for (j = 0; j < i; j++)
            {
                if (gr_poly_equal(f, fac->entries + j, ctx) != T_FALSE)
                {
                    flint_printf("FAIL (repeated factor)\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("algorithm = %d\n", algorithm);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("fac = %{gr_poly*}\n\n", fac->entries, fac->length, ctx);
                    flint_printf("exp = %{fmpz*}\n\n", exp->entries, exp->length);
                    flint_abort();
                }
            }
        }

        /* Test the gr_factor interface on the polynomial ring */
        if (n_randint(state, 4) == 0 && P->length > 0)
        {
            gr_ctx_t pctx;
            gr_vec_t fac2;
            fmpz_vec_t exp2;
            gr_poly_t c2;

            gr_ctx_init_gr_poly(pctx, ctx);
            gr_vec_init(fac2, 0, pctx);
            fmpz_vec_init(exp2, 0);
            gr_poly_init(c2, ctx);

            status = gr_factor(c2, fac2, exp2, P, 0, pctx);

            if (status != GR_SUCCESS || fac2->length != fac->length)
            {
                flint_printf("FAIL (gr_factor)\n\n");
                gr_ctx_println(ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_vec_clear(fac2, pctx);
            fmpz_vec_clear(exp2);
            gr_poly_clear(c2, ctx);
            gr_ctx_clear(pctx);
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_vec_clear(fac, ctx);
        fmpz_vec_clear(exp);
        gr_heap_clear(c, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
