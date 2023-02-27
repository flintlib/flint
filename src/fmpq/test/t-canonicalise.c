/*
    Copyright (C) 2011 Fredrik Johansson

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
#include "fmpq.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("canonicalise....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t x;
        fmpz_t mult;

        fmpq_init(x);
        fmpq_randtest(x, state, 200);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: expected fmpq_randtest output to be canonical\n");
            fmpq_print(x);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_init(mult);
        fmpz_randtest_not_zero(mult, state, 200);
        fmpz_add_ui(mult, mult, UWORD(1));

        fmpz_mul(&x->num, &x->num, mult);
        fmpz_mul(&x->den, &x->den, mult);

        if (fmpq_is_canonical(x))
        {
            flint_printf("FAIL: expected fmpq_is_canonical to detect common factor\n");
            fmpq_print(x);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_canonicalise(x);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical after calling fmpq_canonicalise\n");
            fmpq_print(x);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_neg(&x->den, &x->den);

        if (fmpq_is_canonical(x))
        {
            flint_printf("FAIL: negative denominator reported as being canonical\n");
            fmpq_print(x);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_canonicalise(x);

        if (!fmpq_is_canonical(x))
        {
            flint_printf("FAIL: result not canonical after calling fmpq_canonicalise\n");
            fmpq_print(x);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(mult);
        fmpq_clear(x);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
