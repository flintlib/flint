/*
    Copyright (C) 2012 Fredrik Johansson

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

static void
_taylor_simple(fmpz * poly, const fmpz_t c, slong n)
{
    slong i, j;

    for (i = n - 2; i >= 0; i--)
        for (j = i; j < n - 1; j++)
            fmpz_addmul(poly + j, poly + j + 1, c);
}

static void
taylor_simple(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);
    _taylor_simple(g->coeffs, c, g->length);
}

TEST_FUNCTION_START(fmpz_poly_taylor_shift_horner, state)
{
    int i;

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
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n");
            fmpz_poly_print(g); flint_printf("\n");
            fflush(stdout);
            flint_abort();
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
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n");
            fmpz_poly_print(g); flint_printf("\n");
            fmpz_poly_print(h1); flint_printf("\n");
            fmpz_poly_print(h2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h1);
        fmpz_poly_clear(h2);
        fmpz_clear(c);
    }

    /* Shift by 1 and -1 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, h1, h2;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(h1);
        fmpz_poly_init(h2);

        fmpz_init(c);

        if (n_randint(state, 30))
            fmpz_poly_randtest(f, state, 1 + n_randint(state, 100),
                                         1 + n_randint(state, 10000));
        else
            fmpz_poly_randtest(f, state, 1 + n_randint(state, 100),
                                         1 + n_randint(state, 200));

        fmpz_set_si(c, n_randint(state, 2) ? 1 : -1);
        fmpz_poly_taylor_shift_horner(h1, f, c);
        fmpz_neg(c, c);
        taylor_simple(h2, h1, c);

        if (!fmpz_poly_equal(f, h2))
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n");
            fmpz_poly_print(h1); flint_printf("\n");
            fmpz_poly_print(h2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(h1);
        fmpz_poly_clear(h2);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
