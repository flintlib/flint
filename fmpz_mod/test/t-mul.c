/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "profiler.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);
   
    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;   /* p not nec prime */
        fmpz_t a, b, c, d, f1, f2;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(f1);
        fmpz_init(f2);
        fmpz_mod_ctx_init(fpctx, p);

        for (j = 0; j < 20; j++)
        {
            fmpz_randtest_unsigned(p, state, 200);
            fmpz_add_ui(p, p, 1);
            fmpz_mod_ctx_set_mod(fpctx, p);

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            fmpz_randtest_mod(c, state, p);
            fmpz_randtest_mod(d, state, p);
            fmpz_randtest_mod(f1, state, p);
            fmpz_randtest_mod(f2, state, p);

            fmpz_mul(f1, a, b);
            fmpz_mod(f1, f1, p);
            fmpz_mul(f2, b, c);
            fmpz_mod(f2, f2, p);

            fmpz_mod_mul(d, a, b, fpctx);
            if (!fmpz_equal(d, f1))
            {
                printf("FAIL\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
            }

            fmpz_mod_mul(a, a, b, fpctx);
            if (!fmpz_equal(a, f1))
            {
                printf("FAIL\ncheck aliasing first");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
            }

            fmpz_mod_mul(c, b, c, fpctx);
            if (!fmpz_equal(c, f2))
            {
                printf("FAIL\ncheck aliasing second");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
            }

        }

        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(f1);
        fmpz_clear(f2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
