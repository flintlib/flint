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

    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("get/set_coeff_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        fmpq_t x, y;
        long coeff, len;

        fmpq_poly_init(a);
        fmpq_init(x);
        fmpq_init(y);
        len = (long) (n_randint(state, 100) + 1);

        for (j = 0; j < 50; j++)
        {
            fmpq_randtest(x, state, 200);
            coeff = (long) n_randint(state, len);
            fmpq_poly_set_coeff_fmpq(a, coeff, x);
            fmpq_poly_get_coeff_fmpq(y, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (fmpq_equal(x, y) && !cflags);
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("a     = "), fmpq_poly_debug(a), printf("\n\n");
                printf("coeff = %ld\n\n", coeff);
                printf("len   = %ld\n\n", len);
                printf("cflags = %lu\n\n", cflags);
                printf("x = "), fmpq_print(x), printf("\n");
                printf("y = "), fmpq_print(y), printf("\n");
                abort();
            }
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_poly_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
