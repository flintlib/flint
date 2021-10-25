/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("init/clear....");
    fflush(stdout);

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f;
       slong len, coeff_bits, exp_bits;

       fmpz_mpoly_ctx_init_rand(ctx, state, 20);

       fmpz_mpoly_init(f, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
       fmpz_mpoly_assert_canonical(f, ctx);

       fmpz_mpoly_clear(f, ctx);        
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

