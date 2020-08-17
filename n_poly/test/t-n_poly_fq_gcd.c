/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "n_poly.h"


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("n_poly_fq_gcd....");
    fflush(stdout);

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

        for (j = 0; j < 5; j++)
        {
            n_poly_fq_randtest(a, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(b, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(c, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(d, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(e, state, n_randint(state, 15), ctx);
            n_poly_fq_mul(b, b, e, ctx);
            n_poly_fq_mul(c, c, e, ctx);

            n_poly_fq_get_fq_nmod_poly(B, b, ctx);
            n_poly_fq_get_fq_nmod_poly(C, c, ctx);
            fq_nmod_poly_gcd(A, B, C, ctx);
            n_poly_fq_set_fq_nmod_poly(a, A, ctx);

            n_poly_fq_gcd(d, b, c, ctx);
            n_poly_fq_set(e, b, ctx);
            n_poly_fq_gcd(b, b, c, ctx);
            n_poly_fq_gcd(c, e, c, ctx);
            if (!n_poly_fq_equal(a, d, ctx) ||
                !n_poly_fq_equal(a, b, ctx) ||
                !n_poly_fq_equal(a, c, ctx))
            {
                flint_printf("FAIL\n i = %wd, j = %wd\n", i, j);
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


    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
