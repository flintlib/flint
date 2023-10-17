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
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_get_set_coeff_ui, state)
{
    int i, j, result;
    ulong cflags = UWORD(0);

    ulong n;
    fmpq_t n_fmpq;

    fmpq_init(n_fmpq);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        slong coeff, len;

        fmpq_poly_init(a);
        len = (slong) (n_randint(state, 100) + 1);
        fmpq_poly_randtest(a, state, len, 100);

        for (j = 0; j < 1000; j++)
        {
            n = n_randtest(state);
            coeff = n_randint(state, len);
            fmpq_poly_set_coeff_ui(a, coeff, n);
            fmpq_poly_get_coeff_fmpq(n_fmpq, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (fmpz_cmp_ui(fmpq_denref(n_fmpq), 1) == 0
                   && fmpz_cmp_ui(fmpq_numref(n_fmpq), n) == 0
                   && !cflags);
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("a      = "), fmpq_poly_debug(a), flint_printf("\n");
                flint_printf("len    = %wd\n", len);
                flint_printf("coeff  = %wd\n", coeff);
                flint_printf("cflags = %wu\n", cflags);
                flint_printf("n      = %wu\n", n);
                flint_printf("n_fmpq  = "); fmpq_print(n_fmpq); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_poly_clear(a);
    }
    fmpq_clear(n_fmpq);

    TEST_FUNCTION_END(state);
}
