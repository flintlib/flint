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
   
    flint_printf("unity_zp_equal....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        ulong p;
        fmpz_t n;
        unity_zp f, g;

        p = n_randprime(state, 2 + n_randint(state, 6), 0);

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f, p, 1, n);
        unity_zp_init(g, p, 1, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, p);

            fmpz_randtest_not_zero(val, state, 300);
            unity_zp_coeff_set_fmpz(f, ind, val);
            unity_zp_coeff_set_fmpz(g, ind, val);

            fmpz_clear(val);
        }

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

