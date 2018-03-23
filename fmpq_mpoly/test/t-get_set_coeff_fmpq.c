/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_coeff_fmpq....");
    fflush(stdout);

    /* Set coeff and get coeff and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t c, d;
        slong len, coeff_bits, exp_bits, index;

        fmpq_init(c);
        fmpq_init(d);
        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpq_randtest(c, state, n_randint(state, 100) + 1);

            len = f->zpoly->length;
            index = n_randint(state, len + 1);
            fmpq_mpoly_set_coeff_fmpq(f, index, c, ctx);

            if (!fmpq_is_zero(c))
            {
                fmpq_mpoly_get_coeff_fmpq(d, f, index, ctx);
                result = fmpq_equal(c, d);
            } else
            {
                result = (index >= len) || (f->zpoly->length == len - 1);
            }

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Set coeff and get coeff and compare\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_clear(c);
        fmpq_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

