/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_vec_mul, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C;
        TEMPLATE(T, struct) * a, * c;
        TEMPLATE(T, struct) ** aa, ** cc;
        slong j, m, n, alen;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        TEMPLATE(T, mat_init)(C, 1, n, ctx);
        TEMPLATE(T, mat_init)(A, 1, m, ctx);
        TEMPLATE(T, mat_init)(B, m, n, ctx);
        c = _TEMPLATE(T, vec_init)(n, ctx);
        a = _TEMPLATE(T, vec_init)(alen, ctx);

        TEMPLATE(T, mat_randtest)(B, state, ctx);
        _TEMPLATE(T, vec_randtest)(c, state, n, ctx);
        _TEMPLATE(T, vec_randtest)(a, state, alen, ctx);

        cc = FLINT_ARRAY_ALLOC(n, TEMPLATE(T, struct) *);
        for (j = 0; j < n; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, TEMPLATE(T, struct));
            TEMPLATE(T, init)(cc[j], ctx);
            TEMPLATE(T, set)(cc[j], c + j, ctx);
        }

        aa = FLINT_ARRAY_ALLOC(alen, TEMPLATE(T, struct) *);
        for (j = 0; j < alen; j++)
        {
            aa[j] = FLINT_ARRAY_ALLOC(1, TEMPLATE(T, struct));
            TEMPLATE(T, init)(aa[j], ctx);
            TEMPLATE(T, set)(aa[j], a + j, ctx);
        }

        TEMPLATE(T, mat_vec_mul)(c, a, alen, B, ctx);
        TEMPLATE(T, mat_vec_mul_ptr)(cc,
                        (const TEMPLATE(T, struct) * const *)aa, alen, B, ctx);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < m && j < alen; j++)
            TEMPLATE(T, set)(TEMPLATE(T, mat_entry)(A, 0, j), a + j, ctx);

        TEMPLATE(T, mat_mul)(C, A, B, ctx);

        for (j = 0; j < n; j++)
        {
            if (!TEMPLATE(T, equal)(TEMPLATE(T, mat_entry)(C, 0, j), c + j, ctx) ||
                !TEMPLATE(T, equal)(TEMPLATE(T, mat_entry)(C, 0, j), cc[j], ctx))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);
        TEMPLATE(T, mat_clear)(C, ctx);

        _TEMPLATE(T, vec_clear)(c, n, ctx);
        _TEMPLATE(T, vec_clear)(a, alen, ctx);

        for (j = 0; j < n; j++)
        {
            TEMPLATE(T, clear)(cc[j], ctx);
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < alen; j++)
        {
            TEMPLATE(T, clear)(aa[j], ctx);
            flint_free(aa[j]);
        }
        flint_free(aa);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
