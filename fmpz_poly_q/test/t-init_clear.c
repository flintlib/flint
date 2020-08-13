/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("init/clear... ");
    fflush(stdout);

    

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a;
        slong len1 = n_randint(state, 50);
        slong len2 = n_randint(state, 50);
        flint_bitcnt_t bits1 = n_randint(state, 50);
        flint_bitcnt_t bits2 = n_randint(state, 50);

        fmpz_poly_q_init(a);
        fmpz_poly_q_randtest(a, state, len1, bits1, len2, bits2);
        fmpz_poly_q_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
