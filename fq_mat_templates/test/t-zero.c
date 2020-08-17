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

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);

    printf("zero/is_zero....");
    fflush(stdout);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A;
        slong m, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_zero) (A, ctx);

        if (!TEMPLATE(T, mat_is_zero) (A, ctx))
        {
            printf("FAIL: expected matrix to be zero\n");
            abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);
            TEMPLATE(T, randtest_not_zero) (TEMPLATE(T, mat_entry) (A, m, n),
                                            state, ctx);

            if (TEMPLATE(T, mat_is_zero) (A, ctx))
            {
                printf("FAIL: expected matrix not to be zero\n");
                abort();
            }
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    printf("PASS\n");
    return 0;
}


#endif
