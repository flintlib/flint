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
   
    flint_printf("unity_zpq_pow....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        ulong p, q, pow;
        fmpz_t n;
        unity_zpq f, g, h1, h2;

        p = n_randprime(state, 2 + n_randint(state, 6), 0);
        q = n_randprime(state, 2 + n_randint(state, 6), 0);

        pow =  n_randint(state, 32);        

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zpq_init(f, q, p, n);
        unity_zpq_init(g, q, p, n);
        unity_zpq_init(h1, q, p, n);
        unity_zpq_init(h2, q, p, n);

        for (j = 0; j < 100; j++)
        {
            ulong x, y;
            fmpz_t val;

            fmpz_init(val);

            x = n_randint(state, p);
            y = n_randint(state, q);
            
            fmpz_randtest_not_zero(val, state, 200);

            unity_zpq_coeff_set_fmpz(h1, y, x, val);

            fmpz_clear(val);
        }

        unity_zpq_copy(h2, h1);

        unity_zpq_pow_ui(f, h2, pow);
        if (pow == 0)
        {
            unity_zpq_coeff_set_ui(g, 0, 0, 1);
        } else {
            for (j = 0; j < pow; j++)
            {
                unity_zpq_mul(g, h1, h2);
                unity_zpq_swap(h2, g);
            }
        }

        if (unity_zpq_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zpq_clear(f);
        unity_zpq_clear(g);
        unity_zpq_clear(h1);
        unity_zpq_clear(h2);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

