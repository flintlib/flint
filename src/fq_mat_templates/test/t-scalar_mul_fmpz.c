/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_scalar_mul_fmpz, state)
{
    slong rep;
    
    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C, D;
        fmpz_t c;
        slong m, n;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(c);
        TEMPLATE(T, mat_init)(A, m, n, ctx);
        TEMPLATE(T, mat_init)(B, m, n, ctx);
        TEMPLATE(T, mat_init)(C, m, n, ctx);
        TEMPLATE(T, mat_init)(D, m, n, ctx);
        
        fmpz_randtest(c, state, 100);
        TEMPLATE(T, mat_randtest)(A, state, ctx);
        TEMPLATE(T, mat_randtest)(B, state, ctx);

        TEMPLATE(T, mat_scalar_mul_fmpz)(C, A, c, ctx);
        fmpz_sub_ui(c, c, 1);
        TEMPLATE(T, mat_scalar_mul_fmpz)(D, A, c, ctx);

        /* c*A - (c-1)*A == A */
        TEMPLATE(T, mat_sub)(D, C, D, ctx);

        if (!TEMPLATE(T, mat_equal)(A, D, ctx))
            TEST_FUNCTION_FAIL("");

        /* Aliasing */
        TEMPLATE(T, mat_scalar_mul_fmpz)(C, A, c, ctx);
        TEMPLATE(T, mat_scalar_mul_fmpz)(A, A, c, ctx);

        if (!TEMPLATE(T, mat_equal)(A, C, ctx))
            TEST_FUNCTION_FAIL("");

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);
        TEMPLATE(T, mat_clear)(C, ctx);
        TEMPLATE(T, mat_clear)(D, ctx);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
