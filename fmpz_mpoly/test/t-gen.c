/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("gen/is_gen....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f1, f2;
        ordering_t ord;
        slong nvars, len, coeff_bits, exp_bits, k1, k2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f1, ctx);
        fmpz_mpoly_init(f2, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(f1, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(f2, state, len, coeff_bits, exp_bits, ctx);

        k1 = n_randint(state, nvars);
        k2 = n_randint(state, nvars);
        fmpz_mpoly_gen(f1, k1, ctx);
        fmpz_mpoly_assert_canonical(f1, ctx);
        fmpz_mpoly_gen(f2, k2, ctx);
        fmpz_mpoly_assert_canonical(f2, ctx);

        result = 1;
        result = result && fmpz_mpoly_is_gen(f1, k1, ctx);
        result = result && fmpz_mpoly_is_gen(f1, -WORD(1), ctx);
        result = result && fmpz_mpoly_is_gen(f2, k2, ctx);
        result = result && fmpz_mpoly_is_gen(f2, -WORD(1), ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check one generator\ni = %wd\n", i);
            flint_abort();
        }

        fmpz_mpoly_mul_johnson(f1, f1, f2, ctx);
        result = 1;
        result = result && !fmpz_mpoly_is_gen(f1, k1, ctx);
        result = result && !fmpz_mpoly_is_gen(f1, k2, ctx);
        result = result && !fmpz_mpoly_is_gen(f1, -WORD(1), ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check product of two generators\ni = %wd\n", i);
            flint_abort();
        }

        fmpz_mpoly_clear(f1, ctx);
        fmpz_mpoly_clear(f2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

