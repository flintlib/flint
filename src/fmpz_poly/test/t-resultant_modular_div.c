/*
    Copyright (C) 2015 Tommy Hofmann

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

TEST_FUNCTION_START(fmpz_poly_resultant_modular_div, state)
{
    int i, result;

    /* Just one specific test */
    {
        fmpz_poly_t f, g;
        fmpz_t a, b, div;
        slong nbits;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(div);

        fmpz_poly_set_str(f, "11  -15 -2 -2 17 0 0 6 0 -5 1 -1");
        fmpz_poly_set_str(g, "9  2 1 1 1 1 1 0 -1 -2");
        fmpz_set_str(div, "11", 10);
        nbits = 42;
        fmpz_poly_resultant_modular_div(a, f, g, div, nbits);
        /* The result is -44081924855067 = -4007447714097 * 11
         * We supply 11 and the missing divisor is less then 2^35 */
        fmpz_set_str(b, "-4007447714097", 10);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f(x) = "), fmpz_poly_print_pretty(f, "x"), flint_printf("\n\n");
            flint_printf("g(x) = "), fmpz_poly_print_pretty(g, "x"), flint_printf("\n\n");
            flint_printf("res(f, h)/div  = "), fmpz_print(b), flint_printf("\n\n");
            flint_printf("res_mod_div(f, h) = "), fmpz_print(a), flint_printf("\n\n");
            flint_printf("divr = "), fmpz_print(div), flint_printf("\n\n");
            flint_printf("bitsbound = %wd", nbits), flint_printf("\n\n");

            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(div);
    }

    /* Check that R(fg, h) = R(f, h) R(g, h) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;
        fmpz_poly_t f, g, h, p;
        slong nbits;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(p);
        fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 50), 100);

        fmpz_poly_resultant_modular(a, f, h);
        fmpz_poly_resultant_modular(b, g, h);

        if (fmpz_is_zero(b) || fmpz_is_zero(a))
        {
           fmpz_clear(b);
           fmpz_clear(a);
           fmpz_poly_clear(f);
           fmpz_poly_clear(g);
           fmpz_poly_clear(h);
           continue;
        }

        fmpz_mul(c, a, b);
        fmpz_poly_mul(p, f, g);
        nbits = (slong)fmpz_bits(a) + 1; /* for sign */
        fmpz_poly_resultant_modular_div(d, p, h, b, nbits);

        result = (fmpz_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p(x) = "), fmpz_poly_print_pretty(p, "x"), flint_printf("\n\n");
            flint_printf("h(x) = "), fmpz_poly_print_pretty(h, "x"), flint_printf("\n\n");
            flint_printf("res(p, h) = "), fmpz_print(c), flint_printf("\n\n");
            flint_printf("res(p, h) = "), fmpz_print(a), flint_printf(" * "), fmpz_print(b), flint_printf("\n\n");
            flint_printf("supplied divisor = "), fmpz_print(b), flint_printf("\n\n");
            flint_printf("result should be = "), fmpz_print(a), flint_printf("\n\n");
            flint_printf("res(p, h)/div    = "), fmpz_print(d), flint_printf("\n\n");
            flint_printf("bitsbound for result = %wd", nbits), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(p);
    }

    TEST_FUNCTION_END(state);
}
