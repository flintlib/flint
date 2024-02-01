/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_frobenius, state)
{
    slong ix, jx;
    int result;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fq_zech_ctx_t ctx;

        fq_zech_ctx_init_randtest(ctx, state, 2);

        /* Check aliasing: a = frob(a, e) */
        for (jx = 0; jx < 10; jx++)
        {
            fq_zech_t a, b;
            fmpz_t e;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fmpz_init(e);

            fq_zech_randtest(a, state, ctx);
            fmpz_randtest_unsigned(e, state, 6);

            fq_zech_frobenius(b, a, fmpz_get_ui(e), ctx);
            fq_zech_frobenius(a, a, fmpz_get_ui(e), ctx);

            result = (fq_zech_equal(a, b, ctx));
            if (!result)
            {
                flint_printf("FAIL (alias):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fmpz_clear(e);
        }

        /* Compare with exponentiation, for integral values */
        for (jx = 0; jx < 10; jx++)
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

            fq_zech_frobenius(b, a, fmpz_get_ui(e), ctx);

            fmpz_ui_pow_ui(e, fq_zech_ctx_prime(ctx), fmpz_get_ui(e));
            fq_zech_pow(c, a, e, ctx);

            result = (fq_zech_equal(b, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (cmp with pow):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("e = "), fmpz_print(e), flint_printf("\n");
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
