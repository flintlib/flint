/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, assign, state)
{
    int i;

    /* Check that gen does not segfault */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) x;

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 2);
#else
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);
#endif

        TEMPLATE(T, init)(x, ctx);

        TEMPLATE(T, gen)(x, ctx);

        TEMPLATE(T, clear)(x, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
