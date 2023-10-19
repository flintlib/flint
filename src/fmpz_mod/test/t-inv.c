/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"

TEST_FUNCTION_START(fmpz_mod_inv, state)
{
    slong i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_t a, ainv, b, g;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_init(a);
        fmpz_init(ainv);
        fmpz_init(b);
        fmpz_init(g);
        fmpz_mod_ctx_init(fpctx, p);

        for (j = 0; j < 20; j++)
        {
            fmpz_randtest_unsigned(p, state, 200);
            fmpz_add_ui(p, p, 1);
            fmpz_mod_ctx_set_modulus(fpctx, p);

            fmpz_randtest_mod(a, state, p);
            fmpz_gcd(g, a, p);
            if (!fmpz_is_one(g))
            {
                continue;
            }

            fmpz_mod_inv(ainv, a, fpctx);
            fmpz_mod_mul(b, a, ainv, fpctx);
            if (!fmpz_mod_is_one(b, fpctx))
            {
                printf("FAIL\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_inv(a, a, fpctx);
            if (!fmpz_equal(a, ainv))
            {
                printf("FAIL\ncheck aliasing");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(ainv);
        fmpz_clear(b);
        fmpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
