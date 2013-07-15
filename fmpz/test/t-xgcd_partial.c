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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("xgcd_partial....");
    fflush(stdout);

    flint_randinit(state);

    /* Test co2*r1 - co1*r2 = r2_orig */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t co1, co2, f, g, t1, t2, L;

        fmpz_init(co1);
        fmpz_init(co2);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(L);
        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_add_ui(g, g, 1);
        fmpz_randm(f, state, g);
        fmpz_randtest_unsigned(L, state, 200);
        
        fmpz_set(t2, g);
        fmpz_abs(t2, t2);

        fmpz_xgcd_partial(co2, co1, g, f, L);

        fmpz_mul(t1, co2, f);
        fmpz_submul(t1, co1, g);
        fmpz_abs(t1, t1);

        result = fmpz_equal(t1, t2);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("co1 = "), fmpz_print(co1), printf("\n");
            printf("co2 = "), fmpz_print(co2), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("L = "), fmpz_print(L), printf("\n");
            printf("t1 = "), fmpz_print(t1), printf("\n");
            printf("t2 = "), fmpz_print(t2), printf("\n");
            abort();
        }

        fmpz_clear(co1);
        fmpz_clear(co2);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(L);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

