/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"
#include "nfloat.h"

#define NUM_RINGS 17

static int
init_ring(gr_ctx_t ctx, int which)
{
    fmpz_t n;
    int status = GR_SUCCESS;

    fmpz_init(n);

    switch (which)
    {
        case 0: gr_ctx_init_fmpz(ctx); break;
        case 1: gr_ctx_init_fmpq(ctx); break;
        case 2: gr_ctx_init_nmod(ctx, 101); break;
        case 3: status = gr_ctx_init_nmod8(ctx, 7); break;
        case 4: status = gr_ctx_init_nmod32(ctx, 1000003); break;
        case 5:
            fmpz_set_ui(n, 1000003);
            gr_ctx_init_fmpz_mod(ctx, n);
            break;
        case 6:
            fmpz_one(n);
            fmpz_mul_2exp(n, n, 100);
            fmpz_add_ui(n, n, 277);
            status = gr_ctx_init_mpn_mod(ctx, n);
            break;
        case 7:
            fmpz_set_ui(n, 5);
            gr_ctx_init_fq(ctx, n, 2, "a");
            break;
        case 8: gr_ctx_init_fq_nmod(ctx, 5, 2, "a"); break;
        case 9: gr_ctx_init_fq_zech(ctx, 5, 2, "a"); break;
        case 10: gr_ctx_init_real_arb(ctx, 64); break;
        case 11: gr_ctx_init_complex_acb(ctx, 64); break;
        case 12: gr_ctx_init_real_float_arf(ctx, 64); break;
        case 13: gr_ctx_init_complex_float_acf(ctx, 64); break;
        case 14: status = nfloat_ctx_init(ctx, 128, 0); break;
        case 15: status = nfloat_complex_ctx_init(ctx, 128, 0); break;
        default: gr_ctx_init_complex_ca(ctx); break;
    }

    fmpz_clear(n);
    return status;
}

/* the result of multiplying a x b by c x d into e x f */
static int
mul_shapes(gr_ctx_t ctx, slong a, slong b, slong c, slong d, slong e, slong f)
{
    gr_mat_t A, B, C;
    int status;

    gr_mat_init(A, a, b, ctx);
    gr_mat_init(B, c, d, ctx);
    gr_mat_init(C, e, f, ctx);
    status = gr_mat_mul(C, A, B, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

/* gr_mat_mul must report incompatible shapes with GR_DOMAIN, also for
   rings that overload matrix multiplication */
TEST_FUNCTION_START(gr_mat_mul, state)
{
    gr_ctx_t ctx;
    int which;

    for (which = 0; which < NUM_RINGS; which++)
    {
        if (init_ring(ctx, which) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("ring %d: init failed\n", which);

        if (mul_shapes(ctx, 2, 2, 2, 2, 3, 3) != GR_DOMAIN
            || mul_shapes(ctx, 2, 3, 2, 3, 2, 3) != GR_DOMAIN
            || mul_shapes(ctx, 30, 40, 30, 40, 30, 40) != GR_DOMAIN
            || mul_shapes(ctx, 40, 30, 30, 40, 40, 30) != GR_DOMAIN
            || mul_shapes(ctx, 40, 30, 30, 40, 30, 40) != GR_DOMAIN
            || mul_shapes(ctx, 0, 3, 2, 0, 0, 0) != GR_DOMAIN)
        {
            gr_ctx_println(ctx);
            TEST_FUNCTION_FAIL("incompatible shapes not detected\n");
        }

        if (mul_shapes(ctx, 3, 2, 2, 4, 3, 4) != GR_SUCCESS
            || mul_shapes(ctx, 3, 0, 0, 4, 3, 4) != GR_SUCCESS)
        {
            gr_ctx_println(ctx);
            TEST_FUNCTION_FAIL("compatible shapes rejected\n");
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
