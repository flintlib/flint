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
   
    flint_printf("fmpz_mont....");
    fflush(stdout);

    /* test mod add */
    for (i = 0; i < 100; i++)
    {
        fmpz_t f1, f2, g, h, n;

        fmpz_init(f1);
        fmpz_init(f2);

        fmpz_init(n);
        fmpz_init(g);
        fmpz_init(h);

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_randtest_unsigned(h, state, 200);

        fmpz_mod(g, g, n);
        fmpz_mod(h, h, n);

        fmpz_add(f1, g, h);
        fmpz_mod(f1, f1, n);

        mod_add(f2, g, h, n);

        if (fmpz_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(f1);
        fmpz_clear(f2);
        fmpz_clear(g);
        fmpz_clear(h);
        fmpz_clear(n);
    }

    /* test mod sub */
    for (i = 0; i < 100; i++)
    {
        fmpz_t f1, f2, g, h, n;

        fmpz_init(f1);
        fmpz_init(f2);

        fmpz_init(n);
        fmpz_init(g);
        fmpz_init(h);

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_randtest_unsigned(h, state, 200);

        fmpz_mod(g, g, n);
        fmpz_mod(h, h, n);

        fmpz_sub(f1, g, h);
        fmpz_mod(f1, f1, n);

        mod_sub(f2, g, h, n);

        if (fmpz_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(f1);
        fmpz_clear(f2);
        fmpz_clear(g);
        fmpz_clear(h);
        fmpz_clear(n);
    }

    /* test mont mod mul */
    for (i = 0; i < 100; i++)
    {
        ulong r_bits;
        fmpz_t f1, f2, g, h, n, ninv, r, rinv, m, t, d;

        fmpz_init(f1);
        fmpz_init(f2);

        fmpz_init(n);
        fmpz_init(ninv);
        fmpz_init(r);
        fmpz_init(rinv);
        fmpz_init(g);
        fmpz_init(h);
        fmpz_init(m);
        fmpz_init(t);
        fmpz_init_set_ui(d, 1);

        fmpz_randtest_unsigned(n, state, 32);
        while (fmpz_equal_ui(n, 0) != 0 || (fmpz_tstbit(n, 0) == 0))
            fmpz_randtest_unsigned(n, state, 32);

        fmpz_randtest_unsigned(g, state, 32);
        fmpz_randtest_unsigned(h, state, 32);

        r_bits = fmpz_bits(n);   
        fmpz_setbit(r, r_bits);

        fmpz_neg(n, n);
        fmpz_xgcd(d, rinv, ninv, r, n);
        fmpz_mod(ninv, ninv, r);
        fmpz_neg(n, n);

        fmpz_mod(g, g, n);
        fmpz_mod(h, h, n);

        fmpz_mul(f1, g, h);
        fmpz_mod(f1, f1, n);

        fmpz_mul_2exp(g, g, r_bits);
        fmpz_mod(g, g, n);
        fmpz_mul_2exp(h, h, r_bits);
        fmpz_mod(h, h, n);

        mod_mul(f2, g, h, n, ninv, m, t, r_bits);
        mod_mul(f2, f2, d, n, ninv, m, t, r_bits);

        if (fmpz_equal(f1, f2) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(f1);
        fmpz_clear(f2);
        fmpz_clear(g);
        fmpz_clear(h);
        fmpz_clear(n);
        fmpz_clear(ninv);
        fmpz_clear(r);
        fmpz_clear(rinv);
        fmpz_clear(d);
        fmpz_clear(m);
        fmpz_clear(t);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

