/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"

TEST_FUNCTION_START(gr_poly_dispersion, state)
{
    for (slong i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        gr_ctx_t ZZ, ctx;
        gr_poly_t f, g, u, t;
        fmpz_t prev_dfg;
        slong prev_length;

        gr_ctx_init_fmpz(ZZ);

        if (n_randint(state, 4))
            gr_ctx_init_fmpz(ctx);
        else
            gr_ctx_init_random(ctx, state);

        fmpz_init(prev_dfg);
        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(u, ctx);
        gr_poly_init(t, ctx);

        int status = GR_SUCCESS;

        /* Pick polynomials with a nontrivial dispersion structure by
         * incrementally adding new factors; perform various consistency checks
         * on the dispersion sets of intermediate results */

        status |= gr_poly_one(f, ctx);
        status |= gr_poly_one(g, ctx);
        prev_length = 0;

        slong li = 1 + n_randint(state, 5);
        for (slong i = 0; i < li && status == GR_SUCCESS; i++)
        {
            for (slong a = 0; a < 10 && gr_poly_is_zero(u, ctx) == T_TRUE; a++)
                gr_poly_randtest(u, state, 2 + n_randint(state, 3), ctx);
            slong lj = 1 + n_randint(state, 5);
            for (slong j = 0; j < lj && status == GR_SUCCESS; j++)
            {
                fmpz_t duu, dfu, dug, dfg;
                fmpz_init(duu);
                fmpz_init(dfu);
                fmpz_init(dug);
                fmpz_init(dfg);

                gr_vec_t disp;
                gr_vec_init(disp, 0, ctx);

                gr_ptr a = gr_heap_init(ctx);

                if (n_randint(state, 2))
                {
                    status |= gr_set_si(a, n_randint(state, 1ul << 20), ctx);
                    status |= gr_poly_taylor_shift(u, u, a, ctx);
                }
                if (n_randint(state, 2))
                    status |= gr_poly_mul(f, f, u, ctx);
                if (n_randint(state, 2))
                    status |= gr_poly_mul(g, g, u, ctx);

                /* overwritten iff nonempty dispersion set */
                fmpz_set_si(dfu, -1);
                /* overwritten whenever successful */
                status |= gr_vec_append(disp, dfu, ZZ);

                status |= gr_poly_dispersion(duu, NULL, u, u, ctx);
                status |= gr_poly_dispersion(dfu, NULL, f, u, ctx);
                status |= gr_poly_dispersion(dug, NULL, u, g, ctx);
                status |= gr_poly_dispersion(dfg, disp, f, g, ctx);

                if (status == GR_UNABLE) {
                    if (ctx->which_ring == GR_CTX_FMPZ)
                        status = GR_TEST_FAIL;
                    goto cleanup_inner;
                }

                if (status == GR_DOMAIN) {
                    if (!(gr_poly_is_zero(f, ctx) == T_TRUE
                          || gr_poly_is_zero(g, ctx) == T_TRUE))
                        status = GR_TEST_FAIL;
                    goto cleanup_inner;
                }

                if (!fmpz_equal(dfg, prev_dfg) && !fmpz_equal(dfg, duu)
                    && !fmpz_equal(dfg, dfu) && !fmpz_equal(dfg, dug))
                    status = GR_TEST_FAIL;

                if (disp->length < prev_length)
                    status = GR_TEST_FAIL;

                if (disp->length > prev_length && fmpz_equal_si(dfu, -1))
                    status = GR_TEST_FAIL;

                if (ctx->which_ring == GR_CTX_FMPZ)
                {
                    for (slong k = 0; k < disp->length; k++)
                    {
                        fmpz_set(a, gr_vec_entry_ptr(disp, k, ZZ));
                        fmpz_poly_taylor_shift((fmpz_poly_struct *) t,
                                               (fmpz_poly_struct *) f, a);
                        fmpz_poly_gcd((fmpz_poly_struct *) t,
                                      (fmpz_poly_struct *) t,
                                      (fmpz_poly_struct *) g);
                        if (fmpz_poly_degree((fmpz_poly_struct *) t) < 1)
                        {
                            status = GR_TEST_FAIL;
                            break;
                        }
                    }
                }

                prev_length = disp->length;

                fmpz_swap(prev_dfg, dfg);

cleanup_inner:
                gr_vec_clear(disp, ZZ);
                fmpz_clear(dfu);
                fmpz_clear(dug);
                fmpz_clear(dfg);
                gr_heap_clear(a, ctx);
            }
        }

        int aliasing_mode = n_randint(state, 8);
        if (aliasing_mode < 2)
            status |= gr_poly_set(g, f, ctx);
        gr_poly_struct * g0 = aliasing_mode == 0 ? f : g;
        gr_poly_struct * g1 = aliasing_mode == 1 ? f : g;

        if (status != GR_SUCCESS || (f->length > 8 && g->length > 8))
            goto epilogue;

        int want_d = n_randint(state, 16);
        int want_disp = n_randint(state, 16);

        fmpz_t d0, d1;
        gr_vec_t disp0, disp1;
        fmpz_init(d0);
        fmpz_init(d1);
        gr_vec_init(disp0, n_randint(state, 4), ZZ);
        gr_vec_init(disp1, n_randint(state, 4), ZZ);

        status |= gr_poly_dispersion_factor(want_d ? d0 : NULL,
                                            want_disp ? disp0 : NULL,
                                            f, g0, ctx);
        status |= gr_poly_dispersion_resultant(want_d ? d1 : NULL,
                                               want_disp ? disp1 : NULL,
                                               f, g1, ctx);

        if (status != GR_SUCCESS) {
            if (ctx->which_ring == GR_CTX_FMPZ ||
                ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR)
                status = GR_TEST_FAIL;
            goto epilogue;
        }

        if (want_d)
            if (!fmpz_equal(d0, d1))
                status = GR_TEST_FAIL;

        if (want_disp)
            if (disp0->length != disp1->length || _gr_vec_equal(disp0->entries, disp1->entries, disp0->length, ZZ) == T_FALSE)
                status = GR_TEST_FAIL;

        if (want_d && want_disp && disp1->length > 0)
            if (!fmpz_equal(d0, gr_vec_entry_ptr(disp1, disp1->length - 1, ZZ)))
                status = GR_TEST_FAIL;

        gr_vec_clear(disp0, ZZ);
        gr_vec_clear(disp1, ZZ);
        fmpz_clear(d0);
        fmpz_clear(d1);

epilogue:

        if (status == GR_TEST_FAIL)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = %{gr_poly}\n", f, ctx);
            flint_printf("g = %{gr_poly}\n", g, ctx);
            flint_printf("u = %{gr_poly}\n", u, ctx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(prev_dfg);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(u, ctx);
        gr_poly_clear(t, ctx);
        gr_ctx_clear(ZZ);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

