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

TEST_FUNCTION_START(fmpz_poly_2norm_normalised_bits, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        fmpz_poly_t f;
        flint_bitcnt_t b1, b2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_poly_init(f);
        do {
            fmpz_poly_randtest(f, state, n_randint(state, 100) + 1, 200);
        } while (f->length == 0);

        fmpz_poly_2norm(a, f);
        fmpz_abs(b, fmpz_poly_lead(f));
        fmpz_fdiv_q(c, a, b);
        b1 = fmpz_bits(c);

        b2 = _fmpz_poly_2norm_normalised_bits(f->coeffs, f->length);

        result = (b1 == b2 || b1 + 1 == b2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            flint_printf("b1 = %wd, b2 = %wd\n", b1, b2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
