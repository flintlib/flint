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
   
    flint_printf("unity_zp_mul....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        ulong p;
        fmpz_t n;
        unity_zp res, left, right, test;

        fmpz_init(n);

        p = n_randprime(state, 2 + n_randint(state, 6), 0);

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(res, p, n);
        unity_zp_init(test, p, n);
        unity_zp_init(left, p, n);
        unity_zp_init(right, p, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val1, val2;

            fmpz_init(val1);
            fmpz_init(val2);

            ind = n_randint(state, p);
            
            fmpz_randtest_not_zero(val1, state, 200);
            fmpz_randtest_not_zero(val2, state, 200);

            unity_zp_coeff_set_fmpz(left, ind, val1);
            unity_zp_coeff_set_fmpz(right, ind, val2);

            fmpz_add(val1, val1, val2);
            unity_zp_coeff_set_fmpz(test, ind, val1);

            fmpz_clear(val1);
            fmpz_clear(val2);
        }

        unity_zp_mul(res, left, right);
        unity_zp_mul(test, right, left);

        if (unity_zp_equal(res, test) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }
    
        fmpz_clear(n);
        unity_zp_clear(res);
        unity_zp_clear(left);
        unity_zp_clear(right);
        unity_zp_clear(test);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

