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
    Copyright (C) 2011 Fredrik Johansson

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

    printf("revert_series_lagrange_fast....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        long n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        do {
            fmpq_poly_randtest(g, state, n_randint(state, 50), 1+n_randint(state,100));
        } while (fmpq_poly_length(g) < 2 || fmpz_is_zero(g->coeffs + 1));
        fmpq_poly_set_coeff_ui(g, 0, 0);
        n = n_randint(state, 50);

        fmpq_poly_revert_series_lagrange_fast(f, g, n);
        fmpq_poly_revert_series_lagrange_fast(g, g, n);

        result = (fmpq_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL (aliasing):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check f(f^(-1)) = id */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        long n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        do {
            if (n_randint(state, 20) == 0)
                fmpq_poly_randtest(g, state,
                    n_randint(state, 50), 1);
            else
                fmpq_poly_randtest(g, state,
                    n_randint(state, 50), 1+n_randint(state,100));
        } while (fmpq_poly_length(g) < 2 || fmpz_is_zero(g->coeffs + 1));
        fmpq_poly_set_coeff_ui(g, 0, 0);
        n = n_randint(state, 50);

        fmpq_poly_revert_series_lagrange_fast(f, g, n);
        fmpq_poly_compose_series(h, g, f, n);

        result = ((n <= 1 && fmpq_poly_is_zero(h)) ||
            (h->length == 2 && fmpz_is_zero(h->coeffs + 0) &&
                fmpz_is_one(h->coeffs + 1)));
        if (!result)
        {
            printf("FAIL (comparison):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            fmpq_poly_print(h), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
