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
    fmpz_t * t;
    FLINT_TEST_INIT(state);
   
    flint_printf("unity_zp_mul2....");
    fflush(stdout);

    t = (fmpz_t*) flint_malloc(sizeof(fmpz_t) * (SQUARING_SPACE));
    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_init(t[i]);

    /* test multiplication in Z[\zeta_2^2] */
    for (i = 0; i < 100; i++)
    {
        ulong p, k;
        fmpz_t n;
        unity_zp f1, f2, g, h;

        p = 2;
        k = 2;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f1, p, k, n);
        unity_zp_init(f2, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(h, p, k, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val1;
            fmpz_t val2;

            fmpz_init(val1);
            fmpz_init(val2);

            ind = n_randint(state, n_pow(p, k));
            
            fmpz_randtest_unsigned(val1, state, 200);
            fmpz_randtest_unsigned(val2, state, 200);

            unity_zp_coeff_set_fmpz(g, ind, val1);
            unity_zp_coeff_set_fmpz(h, ind, val2);

            fmpz_clear(val1);
            fmpz_clear(val2);
        }

        _unity_zp_reduce_cyclotomic(g);
        _unity_zp_reduce_cyclotomic(h);
        unity_zp_mul4(f1, g, h, t);
        unity_zp_mul(f2, g, h);

        if (unity_zp_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }
    
        fmpz_clear(n);
        unity_zp_clear(f1);
        unity_zp_clear(f2);
        unity_zp_clear(g);
        unity_zp_clear(h);
    }

    /* test multiplication in Z[\zeta_2^3] */
    for (i = 0; i < 100; i++)
    {
        ulong p, k;
        fmpz_t n;
        unity_zp f1, f2, g, h;

        p = 2;
        k = 3;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f1, p, k, n);
        unity_zp_init(f2, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(h, p, k, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val1;
            fmpz_t val2;

            fmpz_init(val1);
            fmpz_init(val2);

            ind = n_randint(state, n_pow(p, k));
            
            fmpz_randtest_unsigned(val1, state, 200);
            fmpz_randtest_unsigned(val2, state, 200);

            unity_zp_coeff_set_fmpz(g, ind, val1);
            unity_zp_coeff_set_fmpz(h, ind, val2);

            fmpz_clear(val1);
            fmpz_clear(val2);
        }

        _unity_zp_reduce_cyclotomic(g);
        _unity_zp_reduce_cyclotomic(h);
        unity_zp_mul8(f1, g, h, t);
        unity_zp_mul(f2, g, h);

        if (unity_zp_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }
    
        fmpz_clear(n);
        unity_zp_clear(f1);
        unity_zp_clear(f2);
        unity_zp_clear(g);
        unity_zp_clear(h);
    }

    /* test multiplication in Z[\zeta_2^4] */
    for (i = 0; i < 100; i++)
    {
        ulong p, k;
        fmpz_t n;
        unity_zp f1, f2, g, h;

        p = 2;
        k = 4;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f1, p, k, n);
        unity_zp_init(f2, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(h, p, k, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val1;
            fmpz_t val2;

            fmpz_init(val1);
            fmpz_init(val2);

            ind = n_randint(state, n_pow(p, k));
            
            fmpz_randtest_unsigned(val1, state, 200);
            fmpz_randtest_unsigned(val2, state, 200);

            unity_zp_coeff_set_fmpz(g, ind, val1);
            unity_zp_coeff_set_fmpz(h, ind, val2);

            fmpz_clear(val1);
            fmpz_clear(val2);
        }

        _unity_zp_reduce_cyclotomic(g);
        _unity_zp_reduce_cyclotomic(h);
        unity_zp_mul16(f1, g, h, t);
        unity_zp_mul(f2, g, h);

        if (unity_zp_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }
    
        fmpz_clear(n);
        unity_zp_clear(f1);
        unity_zp_clear(f2);
        unity_zp_clear(g);
        unity_zp_clear(h);
    }

    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_clear(t[i]);
    flint_free(t);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

