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
#include "nmod_poly.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("taylor_shift_horner....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        mp_limb_t c, mod;

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_poly_randtest(f, state, 1 + n_randint(state, 50));
        c = n_randtest(state) % mod;

        nmod_poly_taylor_shift_horner(g, f, c);
        nmod_poly_taylor_shift_horner(f, f, c);

        if (!nmod_poly_equal(g, f))
        {
            printf("FAIL\n");
            nmod_poly_print(f); printf("\n");
            nmod_poly_print(g); printf("\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Compare with composition */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h1, h2;
        mp_limb_t mod, c;

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);
        nmod_poly_init(h1, mod);
        nmod_poly_init(h2, mod);

        nmod_poly_randtest(f, state, 1 + n_randint(state, 50));
        c = n_randtest(state) % mod;

        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, c);

        nmod_poly_taylor_shift_horner(h1, f, c);
        nmod_poly_compose(h2, f, g);

        if (!nmod_poly_equal(h1, h2))
        {
            printf("FAIL\n");
            nmod_poly_print(f); printf("\n");
            nmod_poly_print(g); printf("\n");
            nmod_poly_print(h1); printf("\n");
            nmod_poly_print(h2); printf("\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h1);
        nmod_poly_clear(h2);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
