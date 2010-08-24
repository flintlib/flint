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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_randstate_t state;

    printf("shift_left/right....");
    fflush(stdout);

    fmpq_poly_randinit(state);

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        long shift = (long) n_randint(100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(100), 200);

        fmpq_poly_shift_left(b, a, shift);
        fmpq_poly_shift_left(a, a, shift);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n");
            fmpq_poly_print(b), printf("\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        long shift;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest_not_zero(a, state, n_randint(100) + 1, 200);

        shift = (long) n_randint(a->length);

        fmpq_poly_shift_right(b, a, shift);
        fmpq_poly_shift_right(a, a, shift);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n");
            fmpq_poly_print(b), printf("\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, c;
        long shift = (long) n_randint(100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(100), 200);

        fmpq_poly_shift_left(b, a, shift);
        fmpq_poly_shift_right(c, b, shift);

        result = (fmpq_poly_equal(c, a));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n");
            fmpq_poly_print(b), printf("\n");
            fmpq_poly_print(c), printf("\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    fmpq_poly_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
