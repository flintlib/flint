/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_discriminant, state)
{
    int i, result;

    /* Check that disc(fg) = disc(f) * disc(g) * R(f, g)^2 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, r;
        fmpz_poly_t f, g, p;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(r);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(p);
        do {
           fmpz_poly_randtest(f, state, n_randint(state, 30), 100);
        } while (f->length < 2);
        do {
           fmpz_poly_randtest(g, state, n_randint(state, 30), 100);
        } while (g->length < 2);

        fmpz_poly_discriminant(a, f);
        fmpz_poly_discriminant(b, g);
        fmpz_mul(c, a, b);
        fmpz_poly_mul(p, f, g);
        fmpz_poly_discriminant(d, p);
        fmpz_poly_resultant(r, f, g);
        fmpz_mul(r, r, r);
        fmpz_mul(c, c, r);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("disc(f)  = "), fmpz_print(a), flint_printf("\n\n");
            flint_printf("disc(g)  = "), fmpz_print(b), flint_printf("\n\n");
            flint_printf("disc(fg) = "), fmpz_print(d), flint_printf("\n\n");
            flint_printf("disc(f)*disc(g)*res(f,g) = "), fmpz_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(r);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(p);
    }

    /* Check that disc(f) = 0 for degree < 1 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpz_poly_t f;

        fmpz_init(a);
        fmpz_poly_init(f);

        fmpz_poly_randtest(f, state, 1, 100);

        fmpz_poly_discriminant(a, f);

        result = (fmpz_is_zero(a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("disc(f)  = "), fmpz_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
