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
    int result;
    slong i, j1, j2;
    FLINT_TEST_INIT(state);

    flint_printf("cmp....");
    fflush(stdout);

    /* check polynomial terms are in order */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, mf, mg;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(mf, ctx);
        fmpq_mpoly_init(mg, ctx);

        len = n_randint(state, 100) + 1;
        exp_bits = n_randint(state, 200) + 2;
        exp_bits = n_randint(state, exp_bits) + 2;
        coeff_bits = n_randint(state, 100) + 1;

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_repack_bits(g, f,
                           f->zpoly->bits + n_randint(state, FLINT_BITS), ctx);
        fmpq_mpoly_assert_canonical(f, ctx);
        fmpq_mpoly_assert_canonical(g, ctx);

        for (j1 = 0; j1 < f->zpoly->length; j1++)
        {
        for (j2 = 0; j2 < g->zpoly->length; j2++)
        {
            fmpq_mpoly_get_term_monomial(mf, f, j1, ctx);
            fmpq_mpoly_get_term_monomial(mg, g, j2, ctx);
            result = fmpq_mpoly_cmp(mf, mg, ctx);
            result =   (result == 0 && j1 == j2)
                    || (result == +1 && j1 < j2)
                    || (result == -1 && j1 > j2);
            if (!result)
            {
                flint_printf("FAIL\n"
                             "check polynomial terms are in order\n"
                             "i = %wd, j1 = %wd, j2 = %wd\n", i, j1, j2);
                flint_abort();
            }
        }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(mf, ctx);
        fmpq_mpoly_clear(mg, ctx);
        fmpq_mpoly_ctx_clear(ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
