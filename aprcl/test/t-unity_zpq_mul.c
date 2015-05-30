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
   
    flint_printf("unity_zpq_mul....");
    fflush(stdout);
    
    for (i = 0; i < 1; i++)
    {
        ulong p, q;
        fmpz_t n;
        unity_zpq res, left, right, test;

        fmpz_init(n);

        p = n_randprime(state, 2 + n_randint(state, 6), 0);
        q = n_randprime(state, 2 + n_randint(state, 6), 0);

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zpq_init(res, q, p, n);
        unity_zpq_init(test, q, p, n);
        unity_zpq_init(left, q, p, n);
        unity_zpq_init(right, q, p, n);

        for (j = 0; j < 100; j++)
        {
            ulong x, y;
            fmpz_t val1, val2;

            fmpz_init(val1);
            fmpz_init(val2);

            x = n_randint(state, p);
            y = n_randint(state, q);
            
            fmpz_randtest_not_zero(val1, state, 200);
            fmpz_randtest_not_zero(val2, state, 200);

            unity_zpq_set(left, y, x, val1);
            unity_zpq_set(right, y, x, val2);

            fmpz_clear(val1);
            fmpz_clear(val2);
        }

        unity_zpq_mul(res, left, right);
        unity_zpq_mul(test, right, left);

        if (unity_zpq_equal(res, test) == 0)
        {
            
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zpq_clear(res);
        unity_zpq_clear(left);
        unity_zpq_clear(right);
        unity_zpq_clear(temp);
    }

    flint_printf("PASS\n");
    return 0;
}

