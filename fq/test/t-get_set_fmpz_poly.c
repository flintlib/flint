/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_fmpz_poly... ");
    fflush(stdout);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t x, y;
        fmpz_poly_t z;

        fq_ctx_randtest(ctx, state);
        fq_init(x, ctx);
        fq_init(y, ctx);
        fmpz_poly_init(z);

        for (j = 0; j < 100; j++)
        {
            fq_rand(x, state, ctx);
            fq_rand(y, state, ctx);
            fq_get_fmpz_poly(z, x, ctx);
            fq_set_fmpz_poly(y, z, ctx);

            if (!fq_equal(y, x, ctx))
            {
                flint_printf("FAIL:\n");
                flint_printf("check get/set match i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_poly_clear(z);
        fq_clear(y, ctx);
        fq_clear(x, ctx);
        fq_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

