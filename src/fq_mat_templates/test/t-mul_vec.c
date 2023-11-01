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

TEST_TEMPLATE_FUNCTION_START(T, mat_mul_vec, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C;
        TEMPLATE(T, struct) * b, * c;
        TEMPLATE(T, struct) ** bb, ** cc;
        slong j, m, n, blen;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        blen = n_randint(state, 50);

        TEMPLATE(T, mat_init)(C, m, 1, ctx);
        TEMPLATE(T, mat_init)(A, m, n, ctx);
        TEMPLATE(T, mat_init)(B, n, 1, ctx);
        c = _TEMPLATE(T, vec_init)(m, ctx);
        b = _TEMPLATE(T, vec_init)(blen, ctx);

        TEMPLATE(T, mat_randtest)(A, state, ctx);
        _TEMPLATE(T, vec_randtest)(c, state, m, ctx);
        _TEMPLATE(T, vec_randtest)(b, state, blen, ctx);

        cc = FLINT_ARRAY_ALLOC(m, TEMPLATE(T, struct) *);
        for (j = 0; j < m; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, TEMPLATE(T, struct));
            TEMPLATE(T, init)(cc[j], ctx);
            TEMPLATE(T, set)(cc[j], c + j, ctx);
        }

        bb = FLINT_ARRAY_ALLOC(blen, TEMPLATE(T, struct) *);
        for (j = 0; j < blen; j++)
        {
            bb[j] = FLINT_ARRAY_ALLOC(1, TEMPLATE(T, struct));
            TEMPLATE(T, init)(bb[j], ctx);
            TEMPLATE(T, set)(bb[j], b + j, ctx);
        }

        TEMPLATE(T, mat_mul_vec)(c, A, b, blen, ctx);
        TEMPLATE(T, mat_mul_vec_ptr)(cc, A,
                           (const TEMPLATE(T, struct) * const *)bb, blen, ctx);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < n && j < blen; j++)
            TEMPLATE(T, set)(TEMPLATE(T, mat_entry)(B, j, 0), b + j, ctx);

        TEMPLATE(T, mat_mul)(C, A, B, ctx);

        for (j = 0; j < m; j++)
        {
            if (!TEMPLATE(T, equal)(TEMPLATE(T, mat_entry)(C, j, 0), c + j, ctx) ||
                !TEMPLATE(T, equal)(TEMPLATE(T, mat_entry)(C, j, 0), cc[j], ctx))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        _TEMPLATE(T, vec_clear)(c, m, ctx);
        _TEMPLATE(T, vec_clear)(b, blen, ctx);

        for (j = 0; j < m; j++)
        {
            TEMPLATE(T, clear)(cc[j], ctx);
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < blen; j++)
        {
            TEMPLATE(T, clear)(bb[j], ctx);
            flint_free(bb[j]);
        }
        flint_free(bb);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
