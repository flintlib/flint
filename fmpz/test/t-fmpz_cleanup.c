/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "flint.h"
#include "ulong_extras.h"

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("fmpz_cleanup....");
    fflush(stdout);

    

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        fmpz *A, *B;

        n = n_randint(state, 100);
        A = _fmpz_vec_init(n);
        B = _fmpz_vec_init(n);

        for (i = 0; i < n; i++)
        {
            fmpz_randtest(A + i, state, 1 + n_randint(state, 1000));
            fmpz_randtest(B + i, state, 1 + n_randint(state, 1000));
        }

        for (i = 0; i < n; i++)
        {
            fmpz_addmul(A + n_randint(state, n),
                A + n_randint(state, n), B + n_randint(state, n));
            fmpz_addmul(B + n_randint(state, n),
                A + n_randint(state, n), B + n_randint(state, n));

            if (n_randint(state, 10) == 0)
            {
                
            }
        }

        _fmpz_vec_clear(A, n);
        _fmpz_vec_clear(B, n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

