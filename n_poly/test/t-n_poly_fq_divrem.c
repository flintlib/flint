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

    flint_printf("n_poly_fq_divrem....");
    fflush(stdout);

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fq_nmod_ctx_t ctx;
        n_poly_t a, b, c, d, e, f, q, r;
        fq_nmod_poly_t A, B, Q, R;

        fq_nmod_ctx_randtest(ctx, state);

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);
        n_poly_init(d);
        n_poly_init(e);
        n_poly_init(f);
        n_poly_init(q);
        n_poly_init(r);
        fq_nmod_poly_init(A, ctx);
        fq_nmod_poly_init(B, ctx);
        fq_nmod_poly_init(Q, ctx);
        fq_nmod_poly_init(R, ctx);

        for (j = 0; j < 10; j++)
        {
            n_poly_fq_randtest(a, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(b, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(c, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(d, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(e, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(f, state, n_randint(state, 20), ctx);

            if (n_poly_is_zero(b))
                n_poly_fq_one(b, ctx);

            n_poly_fq_get_fq_nmod_poly(A, a, ctx);
            n_poly_fq_get_fq_nmod_poly(B, b, ctx);
            fq_nmod_poly_divrem(Q, R, A, B, ctx);
            n_poly_fq_set_fq_nmod_poly(q, Q, ctx);
            n_poly_fq_set_fq_nmod_poly(r, R, ctx);

            n_poly_fq_divrem(c, d, a, b, ctx);
            n_poly_fq_divrem(c, d, a, b, ctx);
            if (!n_poly_fq_equal(c, q, ctx) ||
                !n_poly_fq_equal(d, r, ctx))
            {
flint_printf("\n");
flint_printf("q: "); n_poly_fq_print_pretty(q, "x", ctx); flint_printf("\n");
flint_printf("r: "); n_poly_fq_print_pretty(r, "x", ctx); flint_printf("\n");
flint_printf("\n");
flint_printf("c: "); n_poly_fq_print_pretty(c, "x", ctx); flint_printf("\n");
flint_printf("d: "); n_poly_fq_print_pretty(d, "x", ctx); flint_printf("\n");

                flint_printf("FAIL\n i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        n_poly_clear(d);
        n_poly_clear(e);
        n_poly_clear(f);
        n_poly_clear(q);
        n_poly_clear(r);
        fq_nmod_poly_clear(A, ctx);
        fq_nmod_poly_clear(B, ctx);
        fq_nmod_poly_clear(Q, ctx);
        fq_nmod_poly_clear(R, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
