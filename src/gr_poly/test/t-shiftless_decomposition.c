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
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

static int
gr_poly_is_normal(const gr_poly_t pol, gr_ctx_t ctx)
{
    if (pol->length == 0)
        return 1;
    gr_ptr lc = gr_vec_entry_ptr(pol->coeffs, pol->length - 1, ctx);
    return gr_is_zero(lc, ctx) == T_FALSE;
}

TEST_FUNCTION_START(gr_poly_shiftless_decomposition, state)
{
    for (slong i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        gr_ctx_t ctx;
        gr_poly_t f, u;

        if (n_randint(state, 2))
            gr_ctx_init_fmpz(ctx);
        else
            gr_ctx_init_random(ctx, state);
        int calcium = ctx->methods == _ca_methods;

        gr_ptr a = gr_heap_init(ctx);
        gr_poly_init(f, ctx);
        gr_poly_init(u, ctx);

        int status = GR_SUCCESS;

        status |= gr_poly_one(f, ctx);
        slong li = 1 + n_randint(state, calcium ? 2 : 5);

        for (slong i = 0; i < li && status == GR_SUCCESS; i++)
        {
            for (slong j = 0; j < 10 && gr_poly_is_zero(u, ctx) == T_TRUE; j++)
                gr_poly_randtest(u, state, 2 + (calcium ? 0 : n_randint(state, 3)), ctx);
            slong lj = 1 + n_randint(state, calcium ? 5 : 2);
            for (slong j = 0; j < lj && status == GR_SUCCESS; j++)
            {
                if (n_randint(state, 2))
                {
                    status |= gr_set_si(a, n_randint(state, 1ul << 20), ctx);
                    status |= gr_poly_taylor_shift(u, u, a, ctx);
                }
                if (n_randint(state, 2))
                    status |= gr_poly_mul(f, f, u, ctx);
            }
        }

        gr_ptr c;
        gr_ctx_t Pol, ZZ, ZZvec;
        gr_vec_t slfac, slshifts, slmult;

        gr_ctx_init_fmpz(ZZ);
        gr_ctx_init_vector_gr_vec(ZZvec, ZZ);
        gr_ctx_init_gr_poly(Pol, ctx);

        c = gr_heap_init(ctx);
        gr_vec_init(slfac, n_randint(state, 4), Pol);
        gr_vec_init(slshifts, n_randint(state, 4), ZZvec);
        gr_vec_init(slmult, n_randint(state, 4), ZZvec);

        status |= gr_poly_shiftless_decomposition(c, slfac, slshifts, slmult, f, ctx);

        if (status == GR_SUCCESS)
        {
            gr_poly_t f1, gcd, shifted, sqf;

            gr_poly_init(f1, ctx);
            gr_poly_init(gcd, ctx);
            gr_poly_init(shifted, ctx);
            gr_poly_init(sqf, ctx);

            status |= gr_poly_set_scalar(f1, c, ctx);
            for (slong i = 0; i < slfac->length; i++)
            {
                gr_poly_struct * gi = gr_vec_entry_ptr(slfac, i, Pol);
                gr_vec_struct * shi = gr_vec_entry_ptr(slshifts, i, ZZvec);
                gr_vec_struct * mi = gr_vec_entry_ptr(slmult, i, ZZvec);

                if (shi->length != mi->length)
                    status = GR_TEST_FAIL;

                if ((status | gr_poly_squarefree_part(sqf, gi, ctx)) == GR_SUCCESS
                        && gr_poly_equal(sqf, gi, ctx) == T_FALSE)
                    status = GR_TEST_FAIL;

                for (slong i1 = 0; i1 < i; i1++)
                {
                    gr_poly_struct * gi1 = gr_vec_entry_ptr(slfac, i1, Pol);
                    if ((status | gr_gcd(gcd, gi, gi1, Pol)) == GR_SUCCESS)
                        if (gcd->length > 1 && gr_poly_is_normal(gcd, ctx))
                            status = GR_TEST_FAIL;
                }

                for (slong j = 0; j < shi->length; j++)
                {
                    fmpz * shij = gr_vec_entry_ptr(shi, j, ZZ);
                    fmpz * mij = gr_vec_entry_ptr(mi, j, ZZ);

                    if (j > 0 && fmpz_cmp(gr_vec_entry_ptr(shi, j - 1, ZZ), shij) >= 0)
                        status = GR_TEST_FAIL;

                    if (fmpz_sgn(mij) <= 0)
                        status = GR_TEST_FAIL;

                    status |= gr_set_fmpz(a, shij, ctx);
                    status |= gr_poly_taylor_shift(shifted, gi, a, ctx);

                    if (!fmpz_is_zero(shij))
                        if ((status | gr_gcd(gcd, gi, shifted, Pol)) == GR_SUCCESS)
                            if (gcd->length > 1 && gr_poly_is_normal(gcd, ctx))
                                status = GR_TEST_FAIL;

                    status |= gr_pow_fmpz(shifted, shifted, mij, Pol);
                    status |= gr_poly_mul(f1, f1, shifted, ctx);
                }
            }

            if (status == GR_TEST_FAIL)
            {
                flint_printf("FAIL (properties of factors):\n");
                flint_printf("f = %{gr_poly}\n", f, ctx);
                fflush(stdout);
                flint_abort();
            }

            if (status == GR_SUCCESS && gr_poly_equal(f1, f, ctx) == T_FALSE)
            {
                flint_printf("FAIL (equality):\n");
                flint_printf("f = %{gr_poly}\n", f, ctx);
                flint_printf("f1 = %{gr_poly}\n", f1, ctx);
                fflush(stdout);
                flint_abort();
            }

            gr_poly_clear(f1, ctx);
            gr_poly_clear(gcd, ctx);
            gr_poly_clear(shifted, ctx);
            gr_poly_clear(sqf, ctx);
        }
        else if (status == GR_DOMAIN)
        {
            if (gr_ctx_is_finite_characteristic(ctx) == T_FALSE
                && gr_ctx_is_integral_domain(ctx) == T_TRUE
                && gr_poly_is_zero(f, ctx) == T_FALSE)
            {
                flint_printf("FAIL (domain):\n");
                flint_printf("f = %{gr_poly}\n", f, ctx);
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            if (ctx->which_ring == GR_CTX_FMPZ ||
                ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR)
            {
                flint_printf("FAIL (unexpected failure):\n");
                flint_printf("f = %{gr_poly}\n", f, ctx);
                fflush(stdout);
                flint_abort();
            }
        }

        gr_vec_clear(slmult, ZZvec);
        gr_vec_clear(slshifts, ZZvec);
        gr_vec_clear(slfac, Pol);
        gr_heap_clear(c, ctx);

        gr_heap_clear(a, ctx);
        gr_poly_clear(f, ctx);
        gr_poly_clear(u, ctx);
        gr_ctx_clear(Pol);
        gr_ctx_clear(ctx);
        gr_ctx_clear(ZZvec);
        gr_ctx_clear(ZZ);
    }

    TEST_FUNCTION_END(state);
}
