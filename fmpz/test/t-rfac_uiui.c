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

    Copyright (C) 2012 Fredrik Johansson

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

    printf("rfac_uiui... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check rf(x,a) * rf(x+a,b) = rf(x,a+b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t xa, r1, r2, r1r2, r3;
        mp_limb_t x;
        ulong a, b;

        fmpz_init(xa);
        fmpz_init(r1);
        fmpz_init(r2);
        fmpz_init(r1r2);
        fmpz_init(r3);

        x = n_randlimb(state);
        a = n_randint(state, 100);
        b = n_randint(state, 100);

        fmpz_set_ui(xa, x);
        fmpz_add_ui(xa, xa, a);

        fmpz_rfac_uiui(r1, x, a);
        fmpz_rfac_ui(r2, xa, b);
        fmpz_rfac_uiui(r3, x, a+b);

        fmpz_mul(r1r2, r1, r2);

        result = fmpz_equal(r1r2, r3);

        if (!result)
        {
            printf("FAIL\n\n");
            printf("x: %lu", x); printf("\n\n");
            printf("a = %lu, b = %lu\n\n", a, b);
            printf("rf(x,a): "); fmpz_print(r1); printf("\n\n");
            printf("rf(x+a,b): "); fmpz_print(r2); printf("\n\n");
            printf("rf(x,a+b): "); fmpz_print(r3); printf("\n\n");
            abort();
        }

        fmpz_clear(xa);
        fmpz_clear(r1);
        fmpz_clear(r2);
        fmpz_clear(r1r2);
        fmpz_clear(r3);
    }

    flint_randclear(state);
    printf("PASS\n");
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
