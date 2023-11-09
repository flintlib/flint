/*
    Copyright (C) 2023 Mathieu Gouttenoire
    
    This file is part of FLINT.
    
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpzi.h"

TEST_FUNCTION_START(fmpzi_is_prime, state)
{
    fmpzi_t n;
    
    fmpzi_init(n);
    
    int small_gaussian_primes[10][10] = {
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 0},
        {0, 1, 1, 0, 1, 0, 1, 0, 0, 0},
        {0, 1, 0, 1, 0, 1, 0, 1, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 1, 0, 1, 0, 1, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}
    };
    
    for (int a = 0; a < 10; a ++) {
        for (int b = 0; b < 10; b ++) {
            
            fmpzi_set_si_si(n, a, b);
            
            if (fmpzi_is_prime(n) != small_gaussian_primes[a][b]) {
                flint_printf("FAIL\n");
                flint_printf("n = "); fmpzi_print(n); printf("\n");
                flint_abort();
            }
        }
    }
    
    fmpzi_clear(n);
    
    TEST_FUNCTION_END(state);
}
