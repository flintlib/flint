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
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("primitive_part....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 200);

        fmpz_poly_primitive_part(f, g);
        fmpz_poly_primitive_part(g, g);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(g), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    /* Check that content(f) primitive_part(f) = sgn(lead(f)) f */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(c);
        fmpz_poly_randtest_not_zero(f, state, n_randint(state, 100) + 1, 200);

        fmpz_poly_content(c, f);
        if (fmpz_sgn(f->coeffs + f->length - 1) < 0)
            fmpz_neg(c, c);
        fmpz_poly_primitive_part(g, f);
        fmpz_poly_scalar_mul_fmpz(g, g, c);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(g), printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
