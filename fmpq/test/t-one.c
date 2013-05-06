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
#include "fmpz.h"
#include "fmpq.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("one....");
    fflush(stdout);

    /* x == 1 * x */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmpq_randtest(x, state, 200);
        fmpq_one(y);

        fmpq_mul(z, y, x);

        result = fmpq_is_canonical(z) && fmpq_equal(x, z);
        if (!result)
        {
            printf("FAIL:\n");
            printf("x = "), fmpq_print(x), printf("\n");
            printf("y = "), fmpq_print(y), printf("\n");
            printf("z = "), fmpq_print(z), printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* x/x == 1 */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);

        while (fmpq_is_zero(x))
            fmpq_randtest(x, state, 200);

        fmpq_div(y, x, x);

        result = fmpq_is_canonical(y) && fmpq_is_one(y);
        if (!result)
        {
            printf("FAIL:\n");
            printf("x = "), fmpq_print(x), printf("\n");
            printf("y = "), fmpq_print(y), printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

