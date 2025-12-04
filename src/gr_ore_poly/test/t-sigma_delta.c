/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_ore_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_ore_poly_sigma_delta, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;

        gr_ore_poly_ctx_init_randtest2(cctx, ctx, state);

        gr_ptr a = gr_heap_init(cctx);
        gr_ptr b = gr_heap_init(cctx);
        gr_ptr c = gr_heap_init(cctx);
        gr_ptr sa = gr_heap_init(cctx);
        gr_ptr sb = gr_heap_init(cctx);
        gr_ptr sc = gr_heap_init(cctx);
        gr_ptr sc1 = gr_heap_init(cctx);
        gr_ptr da = gr_heap_init(cctx);
        gr_ptr db = gr_heap_init(cctx);
        gr_ptr dc = gr_heap_init(cctx);
        gr_ptr dc1 = gr_heap_init(cctx);

        int status = GR_SUCCESS;

        status |= gr_randtest(a, state, cctx);
        status |= gr_randtest(b, state, cctx);

        int call = n_randint(state, 4);
        switch (call)
        {
            case 0:
                status |= gr_ore_poly_sigma_delta(sa, da, a, ctx);
                break;
            case 1:
                status |= gr_ore_poly_sigma(sa, a, ctx);
                status |= gr_ore_poly_delta(da, a, ctx);
                break;
            case 2:
                status |= gr_set(sa, a, cctx);
                status |= gr_ore_poly_sigma_delta(sa, da, sa, ctx);
                break;
            case 3:
                status |= gr_set(da, a, cctx);
                status |= gr_ore_poly_sigma_delta(sa, da, da, ctx);
                break;
            default:
                FLINT_ASSERT(0);
        }

        status |= gr_ore_poly_sigma_delta(sb, db, b, ctx);

        if (status != GR_SUCCESS)
        {
            if (cctx->which_ring == GR_CTX_GR_POLY)
            {
                switch (POLYNOMIAL_CTX(cctx)->base_ring->which_ring)
                {
                    case GR_CTX_FMPZ:
                    case GR_CTX_FMPQ:
                        break;
                    default:
                        goto cleanup;
                }
                switch (GR_ORE_POLY_CTX(ctx)->which_algebra)
                {
                    case ORE_ALGEBRA_COMMUTATIVE:
                    case ORE_ALGEBRA_DERIVATIVE:
                    case ORE_ALGEBRA_EULER_DERIVATIVE:
                    case ORE_ALGEBRA_FORWARD_SHIFT:
                    case ORE_ALGEBRA_FORWARD_DIFFERENCE:
                    case ORE_ALGEBRA_BACKWARD_SHIFT:
                    case ORE_ALGEBRA_BACKWARD_DIFFERENCE:
                    case ORE_ALGEBRA_Q_SHIFT:
                    case ORE_ALGEBRA_MAHLER:
                        flint_printf("FAIL: unexpected failure\n");
                        flint_abort();
                    default:
                        goto cleanup;
                }
            }
            goto cleanup;
        }

        status = gr_add(c, a, b, cctx);
        status |= gr_ore_poly_sigma_delta(sc, dc, c, ctx);

        if (status == GR_SUCCESS)
        {
            status = gr_add(sc1, sa, sb, cctx);
            if (status == GR_SUCCESS && gr_equal(sc, sc1, cctx) == T_FALSE)
            {
                flint_printf("FAIL: σ(a + b) = σ(a) + σ(b)\n");
                flint_abort();
            }

            status = gr_add(dc1, da, db, cctx);
            if (status == GR_SUCCESS && gr_equal(dc, dc1, cctx) == T_FALSE)
            {
                flint_printf("FAIL: δ(a + b) = δ(a) + δ(b)\n");
                flint_abort();
            }
        }

        status = gr_mul(c, a, b, cctx);
        status |= gr_ore_poly_sigma_delta(sc, dc, c, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mul(sc1, sa, sb, cctx);
            if (status == GR_SUCCESS && gr_equal(sc, sc1, cctx) == T_FALSE)
            {
                flint_printf("FAIL: σ(a·b) = σ(a)·σ(b)\n");
                flint_abort();
            }

            status |= gr_mul(dc1, da, b, cctx);
            status |= gr_addmul(dc1, sa, db, cctx);
            if (status == GR_SUCCESS && gr_equal(dc, dc1, cctx) == T_FALSE)
            {
                flint_printf("FAIL: δ(a·b) = δ(a)·b + σ(a)·δ(b)\n");
                flint_abort();
            }
        }

cleanup:
        gr_heap_clear(dc1, cctx);
        gr_heap_clear(dc, cctx);
        gr_heap_clear(db, cctx);
        gr_heap_clear(da, cctx);
        gr_heap_clear(sc1, cctx);
        gr_heap_clear(sc, cctx);
        gr_heap_clear(sb, cctx);
        gr_heap_clear(sa, cctx);
        gr_heap_clear(c, cctx);
        gr_heap_clear(b, cctx);
        gr_heap_clear(a, cctx);
        gr_ore_poly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
