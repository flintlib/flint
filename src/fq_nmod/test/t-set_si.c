/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_nmod.h"

/* fq_nmod_set_si against fq_nmod_set_fmpz, including WORD_MIN, WORD_MAX and
   the negative multiples of p (which must give a reduced zero) */
TEST_FUNCTION_START(fq_nmod_set_si, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_ctx_t ctx;
        fq_nmod_t a, b;
        fmpz_t f;
        slong x;
        ulong p;

        fq_nmod_ctx_init_randtest(ctx, state, 0);
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fmpz_init(f);
        p = fq_nmod_ctx_prime(ctx);

        switch (n_randint(state, 6))
        {
            case 0: x = WORD_MIN; break;
            case 1: x = WORD_MAX; break;
            case 2: x = (slong) (-(p * (1 + n_randint(state, 100)))); break;
            case 3: x = (slong) (-n_randint(state, 3 * p)); break;
            default: x = z_randtest(state); break;
        }

        fq_nmod_set_si(a, x, ctx);
        fmpz_set_si(f, x);
        fq_nmod_set_fmpz(b, f, ctx);

        if (!fq_nmod_equal(a, b, ctx) || (a->length > 0 && a->coeffs[0] >= p))
            TEST_FUNCTION_FAIL("x = %wd, p = %wu, a: length %wd, constant coefficient %wu\n",
                               x, p, a->length, a->length > 0 ? a->coeffs[0] : UWORD(0));

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fmpz_clear(f);
        fq_nmod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
