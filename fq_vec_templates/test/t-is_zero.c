/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    printf("is_zero....");
    fflush(stdout);

    /* Check zero vector */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a;
        slong len = n_randint(state, 100);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        _TEMPLATE(T, vec_zero) (a, len, ctx);

        result = (_TEMPLATE(T, vec_is_zero) (a, len, ctx));
        if (!result)
        {
            printf("FAIL1:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check non-zero vector */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a;
        slong len = n_randint(state, 100) + 1;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        TEMPLATE(T, one) (a + (len - 1), ctx);

        result = (!_TEMPLATE(T, vec_is_zero) (a, len, ctx));
        if (!result)
        {
            printf("FAIL2:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    printf("PASS\n");
    return 0;
}


#endif
