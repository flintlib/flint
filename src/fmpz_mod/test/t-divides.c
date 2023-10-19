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

TEST_FUNCTION_START(fmpz_mod_divides, state)
{
    flint_bitcnt_t max_modulus_bits = 200;
    slong i, j;

    {
        fmpz_t p, a, b, c;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 12);
        fmpz_mod_ctx_init(fpctx, p);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_set_ui(b, 6);
        fmpz_set_ui(c, 3);
        if (!fmpz_mod_divides(a, b, c, fpctx)
            || !(    fmpz_equal_ui(a, 2)
                  || fmpz_equal_ui(a, 6)
                  || fmpz_equal_ui(a, 10) ))
        {
            printf("FAIL\n");
            flint_printf("check 6/3 is 2, 6, or 10 mod 12\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_set_ui(b, 5);
        fmpz_set_ui(c, 3);
        if (fmpz_mod_divides(a, b, c, fpctx))
        {
            printf("FAIL\n");
            flint_printf("check 3 does not divide 5 mod 12\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p; /* p not nec. prime */
        fmpz_t a, b, c, u;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_mod_ctx_init(fpctx, p);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(u);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest_unsigned(p, state, max_modulus_bits);
            fmpz_add_ui(p, p, 1);
            fmpz_mod_ctx_set_modulus(fpctx, p);

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            fmpz_randtest_mod(c, state, p);
            if (fmpz_mod_divides(a, b, c, fpctx))
            {
                fmpz_mod_mul(u, a, c, fpctx);
                if (!fmpz_equal(u, b))
                {
                    printf("FAIL\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            fmpz_randtest_mod(c, state, p);
            fmpz_set(a, b);
            if (fmpz_mod_divides(a, a, c, fpctx))
            {
                fmpz_mod_mul(u, a, c, fpctx);
                if (!fmpz_equal(u, b))
                {
                    printf("FAIL\n");
                    flint_printf("aliasing first\ni = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            fmpz_randtest_mod(c, state, p);
            fmpz_set(a, c);
            if (fmpz_mod_divides(a, b, a, fpctx))
            {
                fmpz_mod_mul(u, a, c, fpctx);
                if (!fmpz_equal(u, b))
                {
                    printf("FAIL\n");
                    flint_printf("aliasing second\ni = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_clear(p);
        fmpz_mod_ctx_clear(fpctx);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(u);
    }

    TEST_FUNCTION_END(state);
}
