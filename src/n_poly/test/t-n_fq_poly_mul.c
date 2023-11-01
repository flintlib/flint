/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "n_poly.h"

TEST_FUNCTION_START(n_fq_poly_mul, state)
{
    slong i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_ctx_t ctx;
        n_poly_t a, b, c, d, e;
        fq_nmod_poly_t A, B, C;

        fq_nmod_ctx_randtest(ctx, state);

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);
        n_poly_init(d);
        n_poly_init(e);
        fq_nmod_poly_init(A, ctx);
        fq_nmod_poly_init(B, ctx);
        fq_nmod_poly_init(C, ctx);

        for (j = 0; j < 10; j++)
        {
            n_fq_poly_randtest(a, state, n_randint(state, 20), ctx);
            n_fq_poly_randtest(b, state, n_randint(state, 20), ctx);
            n_fq_poly_randtest(c, state, n_randint(state, 20), ctx);
            n_fq_poly_randtest(d, state, n_randint(state, 20), ctx);

            n_fq_poly_get_fq_nmod_poly(B, b, ctx);
            n_fq_poly_get_fq_nmod_poly(C, c, ctx);
            fq_nmod_poly_mul(A, B, C, ctx);
            n_fq_poly_set_fq_nmod_poly(a, A, ctx);

            n_fq_poly_mul(d, b, c, ctx);
            n_fq_poly_set(e, b, ctx);
            n_fq_poly_mul(b, b, c, ctx);
            n_fq_poly_mul(c, e, c, ctx);

            if (!n_fq_poly_equal(a, d, ctx) ||
                !n_fq_poly_equal(a, b, ctx) ||
                !n_fq_poly_equal(a, c, ctx))
            {
                flint_printf("FAIL\n i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        n_poly_clear(d);
        n_poly_clear(e);
        fq_nmod_poly_clear(A, ctx);
        fq_nmod_poly_clear(B, ctx);
        fq_nmod_poly_clear(C, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
