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
#include "nmod_vec.h"
#include "nmod_poly.h"
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
        nmod_poly_t f, g;
        mp_limb_t m;
        len_t n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        do {
            nmod_poly_randtest(g, state, n_randint(state, 100));
        } while (nmod_poly_get_coeff_ui(g, 1) == 0);
        nmod_poly_set_coeff_ui(g, 0, 0);
        do {
            n = n_randint(state, 100);
        } while (n >= m);


        nmod_poly_revert_series_lagrange_fast(f, g, n);
        nmod_poly_revert_series_lagrange_fast(g, g, n);

        result = (nmod_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL (aliasing):\n");
            nmod_poly_print(f), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Check f(f^(-1)) = id */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t m;
        len_t n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        nmod_poly_init(h, m);
        do {
            nmod_poly_randtest(g, state, n_randint(state, 100));
        } while (nmod_poly_get_coeff_ui(g, 1) == 0);
        nmod_poly_set_coeff_ui(g, 0, 0);
        do {
            n = n_randint(state, 100);
        } while (n >= m);

        nmod_poly_revert_series_lagrange_fast(f, g, n);
        nmod_poly_compose_series(h, g, f, n);

        result = ((n <= 1 && nmod_poly_is_zero(h)) ||
            (h->length == 2 && h->coeffs[0] == 0 && h->coeffs[1] == 1));
        if (!result)
        {
            printf("FAIL (comparison):\n");
            nmod_poly_print(g), printf("\n\n");
            nmod_poly_print(f), printf("\n\n");
            nmod_poly_print(h), printf("\n\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
