/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, get_set_fmpz, state)
{
    slong ix, jx;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b;
        fmpz_t s, t;

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
#else
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);
#endif
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        fmpz_init(s);
        fmpz_init(t);

        for (jx = 0; jx < 100 * flint_test_multiplier(); jx++)
        {
            TEMPLATE(T, randtest)(a, state, ctx);
            TEMPLATE(T, randtest)(b, state, ctx);
            fmpz_randbits(t, state, 100);

            TEMPLATE(T, set_fmpz)(a, t, ctx);
            FLINT_TEST(TEMPLATE(T, get_fmpz)(s, a, ctx));

            fmpz_sub(t, t, s);
#if defined(FQ_NMOD_H) || defined(FQ_ZECH_H)
            FLINT_TEST(fmpz_divisible_ui(t, TEMPLATE(T, ctx_prime(ctx))));
#else
            FLINT_TEST(fmpz_divisible(t, TEMPLATE(T, ctx_prime(ctx))));
#endif

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

    TEST_FUNCTION_END(state);
}
#endif
