/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int
main(void)
{
    slong i;
    fmpz_t p;
    int success;

    FLINT_TEST_INIT(state);

    flint_printf("next_smooth_prime....");
    fflush(stdout);
   
    flint_randinit(state);

    fmpz_init_set_ui(p, 2);

    for (i = 0; i < 2000; i++)
    {
        success = fmpz_next_smooth_prime(p, p);
        if (!success)
        {
            break;
        }
        if (1 != fmpz_is_prime(p))
        {
            printf("FAIL\nprimality test failed p = ");
            fmpz_print(p);
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(p);
    
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
