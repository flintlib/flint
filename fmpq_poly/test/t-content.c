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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("content....");
    fflush(stdout);

    flint_randinit(state);

    /* Check that content(a f) = abs(a) content(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t a, b, c;

        fmpq_poly_init(f);
        fmpq_poly_init(g);

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        fmpq_poly_randtest_not_zero(f, state, n_randint(state, 100) + 1, 100);
        fmpq_randtest_not_zero(a, state, 100);

        fmpq_poly_scalar_mul_fmpq(g, f, a);
        fmpq_poly_content(b, g);
        fmpq_poly_content(c, f);
        fmpq_mul(c, a, c);
        fmpq_abs(c, c);

        result = (fmpq_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            fmpq_print(a), printf("\n\n");
            fmpq_print(b), printf("\n\n");
            fmpq_print(c), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
