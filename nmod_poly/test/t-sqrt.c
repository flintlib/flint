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
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("sqrt... ");
    fflush(stdout);

    flint_randinit(state);

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        int square1, square2;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);

        nmod_poly_randtest(a, state, 1 + n_randint(state, 50));

        if (n_randint(state, 2))
            nmod_poly_mul(a, a, a);

        square1 = nmod_poly_sqrt(b, a);
        square2 = nmod_poly_sqrt(a, a);

        if ((square1 != square2) || (square1 && !nmod_poly_equal(a, b)))
        {
            printf("FAIL: aliasing:\n");
            printf("square1 = %d, square2 = %d\n\n", square1, square2);
            printf("a: "); nmod_poly_print(a); printf("\n\n");
            printf("b: "); nmod_poly_print(b); printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        int square;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(c, mod);

        nmod_poly_randtest(a, state, 1 + n_randint(state, 50));
        nmod_poly_mul(b, a, a);
        square = nmod_poly_sqrt(c, b);

        if (!square)
        {
            printf("FAIL: square reported nonsquare:\n");
            printf("a: "); nmod_poly_print(a); printf("\n\n");
            printf("b: "); nmod_poly_print(b); printf("\n\n");
            printf("c: "); nmod_poly_print(c); printf("\n\n");
            abort();
        }

        nmod_poly_mul(c, c, c);
        if (!nmod_poly_equal(c, b))
        {
            printf("FAIL: sqrt(b)^2 != b:\n");
            printf("a: "); nmod_poly_print(a); printf("\n\n");
            printf("b: "); nmod_poly_print(b); printf("\n\n");
            printf("c: "); nmod_poly_print(c); printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        long j;
        int square;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(c, mod);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 50));
        nmod_poly_mul(b, a, a);

        j = n_randint(state, nmod_poly_length(b));
        b->coeffs[j] = n_randint(state, mod);
        _nmod_poly_normalise(b);

        square = nmod_poly_sqrt(c, b);

        if (square)
        {
            nmod_poly_mul(c, c, c);
            if (!nmod_poly_equal(c, b))
            {
                printf("FAIL: sqrt(b)^2 != b:\n");
                printf("a: "); nmod_poly_print(a); printf("\n\n");
                printf("b: "); nmod_poly_print(b); printf("\n\n");
                printf("c: "); nmod_poly_print(c); printf("\n\n");
                abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
