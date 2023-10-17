/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_poly_xgcd_modular, state)
{
    int i, result;

    /* Check s*f + t*g == r */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 150);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 150);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_mul(s, s, f);
        fmpz_poly_mul(t, t, g);
        fmpz_poly_add(s, s, t);

        result = fmpz_poly_equal_fmpz(s, r);
        if (!result)
        {
           flint_printf("FAIL (check s*f + t*g == r):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("s = "), fmpz_poly_print(s); flint_printf("\n");
           flint_printf("t = "), fmpz_poly_print(t); flint_printf("\n");
           flint_printf("d = "), fmpz_poly_print(d); flint_printf("\n");
           flint_printf("r = "), fmpz_print(r); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check s*f + t*g == r with smaller polynomials */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 30), 50);
            fmpz_poly_randtest(g, state, n_randint(state, 30), 50);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_mul(s, s, f);
        fmpz_poly_mul(t, t, g);
        fmpz_poly_add(s, s, t);

        result = fmpz_poly_equal_fmpz(s, r);
        if (!result)
        {
           flint_printf("FAIL (check small s*f + t*g == r):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("s = "), fmpz_poly_print(s); flint_printf("\n");
           flint_printf("t = "), fmpz_poly_print(t); flint_printf("\n");
           flint_printf("d = "), fmpz_poly_print(d); flint_printf("\n");
           flint_printf("r = "), fmpz_print(r); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check aliasing of s and f */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_xgcd_modular(r, f, t, f, g);

        result = (fmpz_poly_equal(s, f) || fmpz_is_zero(r));
        if (!result)
        {
           flint_printf("FAIL (alias s and f):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("s = "), fmpz_poly_print(s); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check aliasing of s and g */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_xgcd_modular(r, g, t, f, g);

        result = (fmpz_poly_equal(s, g) || fmpz_is_zero(r));
        if (!result)
        {
           flint_printf("FAIL (alias s and g):\n");
           flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("s = "), fmpz_poly_print(s); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check aliasing of t and f */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_xgcd_modular(r, s, f, f, g);

        result = (fmpz_poly_equal(t, f) || fmpz_is_zero(r));
        if (!result)
        {
           flint_printf("FAIL (alias t and f):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("t = "), fmpz_poly_print(t); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check aliasing of t and g */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_t r;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_init(r);

        do {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd_modular(d, f, g);
        } while (d->length != 1);

        fmpz_poly_xgcd_modular(r, s, t, f, g);
        fmpz_poly_xgcd_modular(r, s, g, f, g);

        result = (fmpz_poly_equal(t, g) || fmpz_is_zero(r));
        if (!result)
        {
           flint_printf("FAIL (alias t and g):\n");
           flint_printf("f = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("t = "), fmpz_poly_print(t); flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_clear(r);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
