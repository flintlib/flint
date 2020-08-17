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

#include "templates.h"

#include <stdio.h>
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, rep;
    FLINT_TEST_INIT(state);

    printf("add/sub/neg....");
    fflush(stdout);

    for (rep = 0; rep < 50 * flint_test_multiplier(); rep++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, mat_t) A;
        TEMPLATE(T, mat_t) B;
        TEMPLATE(T, mat_t) C;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);

        TEMPLATE(T, mat_neg) (C, A, ctx);
        TEMPLATE(T, mat_add) (A, A, B, ctx);
        TEMPLATE(T, mat_sub) (A, A, B, ctx);
        TEMPLATE(T, mat_neg) (A, A, ctx);

        if (!TEMPLATE(T, mat_equal) (A, C, ctx))
        {
            printf("FAIL: matrices not equal!\n");
            abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    printf("PASS\n");
    return 0;
}


#endif
