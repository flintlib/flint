/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default.h"

TEST_FUNCTION_START(fq_default_get_set_fmpz, state)
{
    slong i, j;

    for (i = 0; i < 30*flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t a, b;
        fmpz_t s, t, p;
        slong d;

        fmpz_init(p);
        fmpz_randprime(p, state, n_randint(state, 2) ? 60 : 120, 1);
        d = n_randint(state, 4) + 1,
        fq_default_ctx_init(ctx, p, d, "a");

        fq_default_init(a, ctx);
        fq_default_init(b, ctx);
        fmpz_init(s);
        fmpz_init(t);

        for (j = 0; j < 30*flint_test_multiplier(); j++)
        {
            fq_default_randtest(a, state, ctx);
            fq_default_randtest(b, state, ctx);
            fmpz_randbits(t, state, 100);

            fq_default_set_fmpz(a, t, ctx);
            FLINT_TEST(fq_default_get_fmpz(s, a, ctx));

            fmpz_sub(t, t, s);
            FLINT_TEST(fmpz_divisible(t, p));

            if (fq_default_get_fmpz(t, a, ctx))
            {
                fq_default_set_fmpz(b, t, ctx);
                FLINT_TEST(fq_default_equal(a, b, ctx));
            }
        }

        fmpz_clear(s);
        fmpz_clear(t);
        fq_default_clear(a, ctx);
        fq_default_clear(b, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
