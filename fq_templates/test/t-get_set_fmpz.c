/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "test_helpers.h"

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_fmpz... ");
    fflush(stdout);

    for (i = 0; i < 100*flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b;
        fmpz_t s, t;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        fmpz_init(s);
        fmpz_init(t);

        for (j = 0; j < 100*flint_test_multiplier(); j++)
        {
            TEMPLATE(T, randtest)(a, state, ctx);
            TEMPLATE(T, randtest)(b, state, ctx);
            fmpz_randbits(t, state, 100);

            TEMPLATE(T, set_fmpz)(a, t, ctx);
            FLINT_TEST(TEMPLATE(T, get_fmpz)(s, a, ctx));

            fmpz_sub(t, t, s);
            FLINT_TEST(fmpz_divisible(t, TEMPLATE(T, ctx_prime(ctx))));

            if (TEMPLATE(T, get_fmpz(t, a, ctx)))
            {
                TEMPLATE(T, set_fmpz(b, t, ctx));
                FLINT_TEST(TEMPLATE(T, equal)(a, b, ctx));
            }
        }

        fmpz_clear(s);
        fmpz_clear(t);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
