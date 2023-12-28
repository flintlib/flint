/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

TEST_FUNCTION_START(fmpq_poly_resultant_div, state)
{
    int i, result;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, p;
        fmpq_t x, y, z, zz;
        fmpz_t den;
        slong nbits;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(p);
        fmpq_poly_init(h);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_init(zz);

        fmpz_init(den);

        fmpq_poly_randtest(f, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(h, state, n_randint(state, 50), 100);

        fmpz_set(den, fmpq_poly_denref(f));
        fmpq_poly_scalar_mul_fmpz(f, f, den);

        fmpz_set(den, fmpq_poly_denref(g));
        fmpq_poly_scalar_mul_fmpz(g, g, den);

        fmpz_set(den, fmpq_poly_denref(h));
        fmpq_poly_scalar_mul_fmpz(h, h, den);

        fmpq_poly_mul(p, f, g);

        fmpq_poly_resultant(x, f, h);

        if (!fmpz_is_one(fmpq_denref(x)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_resultant(y, g, h);

        if (!fmpz_is_one(fmpq_denref(y)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("h = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("z = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_resultant(z, p, h);

        if (!fmpz_is_one(fmpq_denref(z)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("p = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (fmpq_is_zero(z))
        {
            fmpq_poly_clear(f);
            fmpq_poly_clear(g);
            fmpq_poly_clear(h);
            fmpq_poly_clear(p);

            fmpq_clear(x);
            fmpq_clear(y);
            fmpq_clear(z);
            fmpq_clear(zz);

            fmpz_clear(den);
            continue;
        }

        nbits = (slong)fmpz_bits(fmpq_numref(y)) + 1;

        fmpq_poly_resultant_div(z, p, h, fmpq_numref(x), nbits);
        fmpq_poly_resultant(zz, p, h);

        result = fmpq_equal(z, y);

        if (!result)
        {
            flint_printf("FAIL (res(p, g)/div == res(p, g)/div:\n");
            flint_printf("p = "), fmpq_poly_print_pretty(p, "x"), flint_printf("\n\n");
            flint_printf("h = "), fmpq_poly_print_pretty(h, "x"), flint_printf("\n\n");
            flint_printf("res(p, h) = "), fmpq_print(zz), flint_printf("\n\n");
            flint_printf("res(p, h) = "), fmpq_print(x), flint_printf(" * "), fmpq_print(y), flint_printf("\n\n");
            flint_printf("supplied divisor = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("nbits = %wu\n\n", nbits);
            flint_printf("divisor found = "), fmpq_print(z), flint_printf("\n\n");
            flint_printf("correct result = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(p);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_clear(zz);

        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
#ifdef __GNUC__
# pragma GCC diagnostic pop
#endif
