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

    Copyright (C) 2011 Sebastian Pancratz

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
    
    printf("init_set....");
    fflush(stdout);
    
    flint_randinit(state);

    /* Small integers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_randtest(a, state, FLINT_BITS - 2);
        fmpz_init_set(b, a);

        result = fmpz_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    /* Large integers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_init_set(b, a);

        result = fmpz_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
