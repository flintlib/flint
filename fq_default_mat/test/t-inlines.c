/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include "fq_default_mat.h"
#include "test_helpers.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("inlines....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg;
        fq_default_ctx_t ctx;
        fq_default_t t;
        fq_default_mat_t a, b, c, d, e, f, g;
        slong n = n_randint(state, 4) + 1;

        fmpz_init(p);
        fmpz_randprime(p, state, n_randint(state, 2) ? 60 : 120, 1);
        deg = n_randint(state, 4) + 1;
        fq_default_ctx_init(ctx, p, deg, "a");

        fq_default_init(t, ctx);
        fq_default_mat_init(a, n, n, ctx);
        fq_default_mat_init(b, n, n, ctx);
        fq_default_mat_init(c, n, n, ctx);
        fq_default_mat_init(d, n, n, ctx);
        fq_default_mat_init(e, n, n, ctx);
        fq_default_mat_init(f, n, n, ctx);
        fq_default_mat_init(g, n, n, ctx);

        fq_default_mat_randtest(a, state, ctx);
        fq_default_mat_randtest(b, state, ctx);
        fq_default_mat_randtest(c, state, ctx);
        fq_default_mat_randtest(d, state, ctx);
        fq_default_mat_randtest(e, state, ctx);
        fq_default_mat_randtest(f, state, ctx);
        fq_default_mat_randtest(g, state, ctx);

        fq_default_mat_add(e, a, b, ctx);
        fq_default_mat_mul(d, e, c, ctx);
        fq_default_mat_mul(f, a, c, ctx);
        fq_default_mat_mul(g, b, c, ctx);
        fq_default_mat_add(e, f, g, ctx);
        FLINT_TEST(fq_default_mat_equal(d, e, ctx));

        fq_default_mat_sub(e, a, b, ctx);
        fq_default_mat_mul(d, e, c, ctx);
        fq_default_mat_mul(f, a, c, ctx);
        fq_default_mat_mul(g, b, c, ctx);
        fq_default_mat_sub(e, f, g, ctx);
        FLINT_TEST(fq_default_mat_equal(d, e, ctx));

        fq_default_clear(t, ctx);
        fq_default_mat_clear(a, ctx);
        fq_default_mat_clear(b, ctx);
        fq_default_mat_clear(c, ctx);
        fq_default_mat_clear(d, ctx);
        fq_default_mat_clear(e, ctx);
        fq_default_mat_clear(f, ctx);
        fq_default_mat_clear(g, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }
    
    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
        return 0;
}
