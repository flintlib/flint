/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "ulong_extras.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, ker;
        slong m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            TEMPLATE(T, ctx_randtest) (ctx, state);
            d = n_randint(state, 2 * m * n + 1);

            TEMPLATE(T, mat_init) (A, m, n, ctx);
            TEMPLATE(T, mat_init) (ker, n, n, ctx);
            TEMPLATE(T, mat_init) (B, m, n, ctx);

            TEMPLATE(T, mat_randrank) (A, state, r, ctx);
            /* Densify */
            if (n_randlimb(state) % 2)
                TEMPLATE(T, mat_randops) (A, d, state, ctx);

            nullity = TEMPLATE(T, mat_nullspace) (ker, A, ctx);
            nulrank = TEMPLATE(T, mat_rank) (ker, ctx);

            if (nullity != nulrank)
            {
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                TEMPLATE(T, mat_print_pretty) (A, ctx);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                TEMPLATE(T, mat_print_pretty) (A, ctx);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, mat_mul) (B, A, ker, ctx);

            if (TEMPLATE(T, mat_rank) (B, ctx) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                TEMPLATE(T, mat_print_pretty) (A, ctx);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, mat_clear) (A, ctx);
            TEMPLATE(T, mat_clear) (ker, ctx);
            TEMPLATE(T, mat_clear) (B, ctx);

            TEMPLATE(T, ctx_clear) (ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
#endif
