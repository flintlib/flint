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

TEST_FUNCTION_START(fmpz_mod_add_sub_neg, state)
{
    slong i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;   /* p not nec prime */
        fmpz_t a, b, c, d, e;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);
        fmpz_mod_ctx_init(fpctx, p);

        for (j = 0; j < 10; j++)
        {
            fmpz_randtest_unsigned(p, state, 300);
            fmpz_add_ui(p, p, 1);
            fmpz_mod_ctx_set_modulus(fpctx, p);

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            fmpz_randtest_mod(c, state, p);
            fmpz_randtest_mod(d, state, p);
            fmpz_randtest_mod(e, state, p);

            fmpz_mod_add(c, a, b, fpctx);
            fmpz_mod_assert_canonical(c, fpctx);
            fmpz_mod_sub(c, c, b, fpctx);
            fmpz_mod_assert_canonical(c, fpctx);
            if (!fmpz_equal(c, a))
            {
                printf("FAIL1\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_sub(d, a, b, fpctx);
            fmpz_mod_assert_canonical(d, fpctx);
            fmpz_mod_sub(e, b, a, fpctx);
            fmpz_mod_assert_canonical(e, fpctx);
            fmpz_mod_neg(d, d, fpctx);
            if (!fmpz_equal(d, e))
            {
                printf("FAIL2\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
