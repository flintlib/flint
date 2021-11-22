/*
    Copyright (C) 2015 Vladimir Glazachev

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
#include "aprcl.h"

int main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);
   
    flint_printf("f_table....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong len, q, p, g;
        mp_ptr table;

        len = n_randint(state, 16);
        while (len < 2)
            len = n_randint(state, 16);

        q = n_randprime(state, len, 1);
        g = n_primitive_root_prime(q);
        p = q - 2;
        table = aprcl_f_table(q);

        for (j = 1; j <= p; j++)
        {
            ulong g_powx, g_powfx;
            g_powx = n_powmod(g, j, q);
            g_powfx = n_powmod(g, table[j], q);
            
            if (n_submod(1, g_powx, q) != g_powfx)
            {
                flint_printf("FAIL:\n");
                flint_printf("1 - %wu != %wu mod %wu\n", g_powx, g_powfx, q);
                fflush(stdout);
                flint_abort();
            }
        }

        _nmod_vec_clear(table);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

