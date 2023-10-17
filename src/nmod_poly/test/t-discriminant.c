/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_discriminant, state)
{
    int i, result;

    /* Check disc(fg) == disc(f) * disc(g) * res(f, g)^2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t x, y, z, r;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        nmod_poly_init(h, n);

        do {
           nmod_poly_randtest(f, state, n_randint(state, 200));
        } while (f->length < 2);
        do {
           nmod_poly_randtest(g, state, n_randint(state, 200));
        } while (g->length < 2);

        y = nmod_poly_discriminant(f);
        z = nmod_poly_discriminant(g);
        y = nmod_mul(y, z, f->mod);
        r = nmod_poly_resultant(f, g);
        r = nmod_mul(r, r, f->mod);
        y = nmod_mul(y, r, f->mod);
        nmod_poly_mul(h, f, g);
        x = nmod_poly_discriminant(h);

        result = (x == y);
        if (!result)
        {
            flint_printf("FAIL (disc(fg) == res(f, g)^2 * disc(f) * disc(g):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            nmod_poly_print(h), flint_printf("\n\n");
            flint_printf("x = %wu\n", x);
            flint_printf("y = %wu\n", y);
            flint_printf("z = %wd\n", z);
            flint_printf("n = %wu\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    /* Check disc(f) == 0 for length < 2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f;
        mp_limb_t y;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);

        nmod_poly_randtest(f, state, 1);

        y = nmod_poly_discriminant(f);

        result = (y == 0);
        if (!result)
        {
            flint_printf("FAIL disc(f) == 0 for len f < 2:\n");
            nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("y = %wu\n", y);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
