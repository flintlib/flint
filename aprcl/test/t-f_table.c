/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

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

    for (i = 0; i < 100; i++)
    {
        ulong len, q, p, g;
        mp_ptr table;

        len = n_randint(state, 16);
        while (len < 2)
            len = n_randint(state, 16);

        q = n_randprime(state, len, 1);
        g = n_primitive_root_prime(q);
        p = q - 2;
        table = f_table(q);

        for (j = 1; j <= p; j++)
        {
            ulong g_powx, g_powfx;
            g_powx = n_powmod(g, j, q);
            g_powfx = n_powmod(g, table[j], q);
            
            if (n_submod(1, g_powx, q) != g_powfx)
            {
                flint_printf("FAIL:\n");
                flint_printf("1 - %wu != %wu mod %wu\n", g_powx, g_powfx, q);
                abort();
            }
        }

        _nmod_vec_clear(table);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

