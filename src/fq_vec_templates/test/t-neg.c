/*
    Copyright (C) 2009, 2010 William Hart
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

TEST_TEMPLATE_FUNCTION_START(T, vec_neg, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a, *b;
        slong len = n_randint(state, 100);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        b = _TEMPLATE(T, vec_init) (len, ctx);
        _TEMPLATE(T, vec_randtest) (a, state, len, ctx);

        _TEMPLATE(T, vec_neg) (b, a, len, ctx);
        _TEMPLATE(T, vec_neg) (a, a, len, ctx);

        result = (_TEMPLATE(T, vec_equal) (a, b, len, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            _TEMPLATE(T, vec_print) (b, len, ctx), printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);
        _TEMPLATE(T, vec_clear) (b, len, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check -(-a) == a */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a, *b;
        slong len = n_randint(state, 100);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        b = _TEMPLATE(T, vec_init) (len, ctx);
        _TEMPLATE(T, vec_randtest) (a, state, len, ctx);

        _TEMPLATE(T, vec_neg) (b, a, len, ctx);
        _TEMPLATE(T, vec_neg) (b, b, len, ctx);

        result = (_TEMPLATE(T, vec_equal) (a, b, len, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            _TEMPLATE(T, vec_print) (b, len, ctx), printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);
        _TEMPLATE(T, vec_clear) (b, len, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
