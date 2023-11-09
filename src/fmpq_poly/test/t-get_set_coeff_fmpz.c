/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_get_set_coeff_fmpz, state)
{
    int i, j, result;
    ulong cflags = UWORD(0);

    fmpq_t n1, n2;

    fmpq_init(n1);
    fmpq_init(n2);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        fmpz_t x1, x2;
        slong coeff, len;

        fmpq_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        len = (slong) (n_randint(state, 100) + 1);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest(x1, state, 200);
            fmpz_set(fmpq_numref(n1), x1);
            fmpz_set_si(fmpq_denref(n1), 1);
            coeff = (slong) n_randint(state, len);
            fmpq_poly_set_coeff_fmpz(a, coeff, x1);
            fmpq_poly_get_coeff_fmpq(n2, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (fmpq_equal(n1, n2) && !cflags);
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a     = "), fmpq_poly_debug(a), flint_printf("\n\n");
                flint_printf("coeff = %wd\n\n", coeff);
                flint_printf("len   = %wd\n\n", len);
                flint_printf("cflags = %wu\n\n", cflags);
                flint_printf("n1 = "); fmpq_print(n1); flint_printf("\n\n");
                flint_printf("n2 = "); fmpq_print(n2); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpq_poly_clear(a);
    }

    fmpq_clear(n1);
    fmpq_clear(n2);

    TEST_FUNCTION_END(state);
}
