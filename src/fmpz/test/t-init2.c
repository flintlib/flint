/*
    Copyright (C) 2009 William Hart

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("init/init2/clear....");
    fflush(stdout);
    
    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;

        fmpz_init2(a, n_randint(state, 100));
        fmpz_clear(a);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;

        fmpz_init(a);
        fmpz_randtest(a, state, SMALL_FMPZ_BITCOUNT_MAX);

        _fmpz_promote_val(a);
        _fmpz_demote_val(a);

        fmpz_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
