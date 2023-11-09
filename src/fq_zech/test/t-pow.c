/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_pow, state)
{
    int j, i, result;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = a^e */
        for (i = 0; i < 100; i++)
        {
            fq_zech_t a, b;
            fmpz_t e;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fmpz_init(e);

            fq_zech_randtest(a, state, ctx);
            fmpz_randtest_unsigned(e, state, 6);

            fq_zech_pow(b, a, e, ctx);
            fq_zech_pow(a, a, e, ctx);

            result = (fq_zech_equal(a, b, ctx));
            if (!result)
            {
                flint_printf("FAIL (alias):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fmpz_clear(e);
        }

        /* Compare with multiplication, for integral values */
        for (i = 0; i < 100; i++)
        {
            fq_zech_t a, b, c;
            fmpz_t e, f;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);
            fmpz_init(f);
            fmpz_init(e);

            fq_zech_randtest(a, state, ctx);
            fmpz_randtest_unsigned(e, state, 6);

            fq_zech_pow(b, a, e, ctx);
            fq_zech_one(c, ctx);
            for (fmpz_one(f); fmpz_cmp(f, e) <= 0; fmpz_add_ui(f, f, 1))
            {
                fq_zech_mul(c, c, a, ctx);
            }

            result = (fq_zech_equal(b, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (cmp with mul):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("e = "), fmpz_print(e), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c, ctx);
            fmpz_clear(e);
            fmpz_clear(f);
        }

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
