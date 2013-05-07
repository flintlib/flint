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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("taylor_shift_horner....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(c);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 20),
                                     1 + n_randint(state, 200));

        fmpz_randtest(c, state, n_randint(state, 200));

        fmpz_poly_taylor_shift_horner(g, f, c);
        fmpz_poly_taylor_shift_horner(f, f, c);

        if (!fmpz_poly_equal(g, f))
        {
            printf("FAIL\n");
            fmpz_poly_print(f); printf("\n");
            fmpz_poly_print(g); printf("\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_clear(c);
    }

    /* Compare with composition */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h1, h2;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h1);
        fmpz_poly_init(h2);

        fmpz_init(c);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 20),
                                     1 + n_randint(state, 200));

        fmpz_randtest(c, state, n_randint(state, 200));

        fmpz_poly_set_coeff_ui(g, 1, 1);
        fmpz_poly_set_coeff_fmpz(g, 0, c);

        fmpz_poly_taylor_shift_horner(h1, f, c);
        fmpz_poly_compose(h2, f, g);

        if (!fmpz_poly_equal(h1, h2))
        {
            printf("FAIL\n");
            fmpz_poly_print(f); printf("\n");
            fmpz_poly_print(g); printf("\n");
            fmpz_poly_print(h1); printf("\n");
            fmpz_poly_print(h2); printf("\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h1);
        fmpz_poly_clear(h2);
        fmpz_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
