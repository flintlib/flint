/*
    Copyright (C) 2014 William Hart
    Copyright (C) 2025 Kacper Proniewski

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_poly.h"
#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpq_poly_discriminant, state)
{
    int i, result;

    /* Check that disc(fg) = disc(f) * disc(g) * R(f, g)^2 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++) {
        fmpq_t a, b, c, d, r;
        fmpq_poly_t f, g, p;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);
        fmpq_init(d);
        fmpq_init(r);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(p);
        do {
            fmpq_poly_randtest(f, state, n_randint(state, 30), 100);
        } while (f->length < 2);
        do {
            fmpq_poly_randtest(g, state, n_randint(state, 30), 100);
        } while (g->length < 2);

        fmpq_poly_discriminant(a, f);
        fmpq_poly_discriminant(b, g);
        fmpq_mul(c, a, b);
        fmpq_poly_mul(p, f, g);
        fmpq_poly_discriminant(d, p);
        fmpq_poly_resultant(r, f, g);
        fmpq_mul(r, r, r);
        fmpq_mul(c, c, r);

        result = (fmpq_equal(c, d));
        if (!result) {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("disc(f)  = "), fmpq_print(a), flint_printf("\n\n");
            flint_printf("disc(g)  = "), fmpq_print(b), flint_printf("\n\n");
            flint_printf("disc(fg) = "), fmpq_print(d), flint_printf("\n\n");
            flint_printf("disc(f)*disc(g)*res(f,g) = "), fmpq_print(c),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);
        fmpq_clear(d);
        fmpq_clear(r);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(p);
    }

    /* Check that disc(f) = 0 for degree < 1 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++) {
        fmpq_t a;
        fmpq_poly_t f;

        fmpq_init(a);
        fmpq_poly_init(f);

        fmpq_poly_randtest(f, state, 1, 100);

        fmpq_poly_discriminant(a, f);

        result = (fmpq_is_zero(a));
        if (!result) {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("disc(f)  = "), fmpq_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
