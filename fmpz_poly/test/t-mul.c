/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    int tmul = 200;
#ifdef _WIN32
    tmul = 1;
#endif
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 500);

        fmpz_poly_mul(a, b, c);
        fmpz_poly_mul(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 500);

        fmpz_poly_mul(a, b, c);
        fmpz_poly_mul(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpz_poly_set(c, b);

        fmpz_poly_mul(a, b, c);
        fmpz_poly_mul(d, b, b);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a1, a2, b, c, d;

        fmpz_poly_init(a1);
        fmpz_poly_init(a2);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 100), 500);
        fmpz_poly_randtest(c, state, n_randint(state, 100), 500);
        fmpz_poly_randtest(d, state, n_randint(state, 100), 500);

        fmpz_poly_mul(a1, b, c);
        fmpz_poly_mul(a2, b, d);
        fmpz_poly_add(a1, a1, a2);

        fmpz_poly_add(c, c, d);
        fmpz_poly_mul(a2, b, c);

        result = (fmpz_poly_equal(a1, a2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a1), flint_printf("\n\n");
            fmpz_poly_print(a2), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a1);
        fmpz_poly_clear(a2);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* Check _fmpz_poly_mul directly */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong len1, len2;
        fmpz_poly_t a, b, out1, out2;

        len1 = n_randint(state, 100) + 1;
        len2 = n_randint(state, 100) + 1;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(out1);
        fmpz_poly_init(out2);
        fmpz_poly_randtest(a, state, len1, 200);
        fmpz_poly_randtest(b, state, len2, 200);

        fmpz_poly_mul(out1, a, b);
        fmpz_poly_fit_length(a, a->alloc + n_randint(state, 10));
        fmpz_poly_fit_length(b, b->alloc + n_randint(state, 10));
        a->length = a->alloc;
        b->length = b->alloc;
        fmpz_poly_fit_length(out2, a->length + b->length - 1);
        if (a->length >= b->length)
            _fmpz_poly_mul(out2->coeffs, a->coeffs, a->length,
                                         b->coeffs, b->length);
        else
            _fmpz_poly_mul(out2->coeffs, b->coeffs, b->length,
                                         a->coeffs, a->length);
        _fmpz_poly_set_length(out2, a->length + b->length - 1);
        _fmpz_poly_normalise(out2);

        result = (fmpz_poly_equal(out1, out2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(out1), flint_printf("\n\n");
            fmpz_poly_print(out2), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(out1);
        fmpz_poly_clear(out2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
